import numpy as np
from scipy.stats import norm
from scipy.special import erfc  # at top

def ground_disp_Mair_1993(x, z0, K, D, vl):
    """
    Calculate the 2D Gaussian based ground deformation. 
    The vertical displacement comes from Mair, R. J., Taylor, R. N., & Bracegirdle,
    A. (1993). Subsurface settlement profiles above tunnels in clays. Geotechnique,
    43(2), 315-320.
    The horizontal displacement comes from O'reilly, M. P., & New, B. M. (1982).
    Settlements above tunnels in the United Kingdom-their magnitude and prediction
    (No. Monograph).

    Parameters
    ----------
    x : float or np.array(dtype=float)
        Distance from tunnel center. Unit: m
    z0 : float or np.array(dtype=float)
        Tunnel center depth. Unit: m
    K : float or np.array(dtype=float)
        Trough width parameter. Unit: Unitless
    D : float
        Tunnel diameter. Unit: m
    vl : float or np.array(dtype=float)
        Tunnel volume loss. Unit: the actual quantity, not percentages

    Returns
    -------
    tuple
        A tuple containing:
        - The vertical displacement (float  or np.array(dtype=float)) in the unit
          of m and upward as positive value.
        - The horizontal displacement (float  or np.array(dtype=float)) in the unit
          of m and positive direction is the same with x.
    """
    # Eq. 3 of Mair
    i = K*z0
    # Eq. 6a of Mair
    s_max = 0.313 * vl * D**2 / i
    # EQ. 1 of Mair
    dv = - s_max * np.exp(-np.power(x, 2)/2/np.power(i,2))
    # Eq. 6 of O'Reilly
    dh = dv * x / z0
    # Return the vertical and horizontal disp
    return dv, dh


def ground_disp_Zhao_2023(x, y, z, vl, d, z0, ys, yf, k, delta):
    """
    Calculate the 3D Gaussian based ground deformation used in Zhao and DeJong (2023). 
    This is a combination of Peck 1969, Attewell 1982, Mair 1993 and Camos et al. 2014,2015,2016.

    Parameters
    ----------
    x : float or np.array(dtype=float)
        Coordinate in the x direction. Unit: m
    y : float or np.array(dtype=float)
        Coordinate in the y direction. Unit: m
    z : float or np.array(dtype=float)
        Coordinate in the z direction. Unit: m
    vl : float or np.array(dtype=float)
        Tunnel volume loss. Unit: the actual quantity, not percentages
    d : float
        Tunnel diameter. Unit: m
    z0 : float
        Tunnel center depth. Unit: m
    ys : float
        Horizontal distance from the tunnel face to the origin.
        Unit: m
    yf : float
        Horizontal distance from the tunnel start to the origin.
        Unit: m
    k : float
        Trough width parameter. Unit: Unitless
    delta : float
        Vertical displacement factor. Unit: Unitless

    Returns
    -------
    tuple
        A tuple containing:
        - The vertical displacement (float or np.array(dtype=float)) in the unit of m and upward as positive value.
        - The horizontal displacement in the x direction (float or np.array(dtype=float)) in the unit of m and positive direction is the same with x.
        - The horizontal displacement in the y direction (float or np.array(dtype=float)) in the unit of m and positive direction is the same with y.
    """

    # Calculate max vertical displacement.
    uz_max = - (np.pi * vl * d**2) / (4 * np.sqrt(2 * np.pi) * k * (z0 - np.min(z)))

    y0 = -norm.ppf(delta) * k * z0

    u_z = uz_max * np.exp(- x**2 / (2 * k**2 * (z0 - z)**2)) * (
        norm.cdf((y - (ys + y0)) / (k * (z0 - z))) - 
        norm.cdf((y - yf) / (k * (z0 - z)))
    )
        
    u_x = (x / (z0 - z)) * u_z

    u_y = vl * d**2 / (8 * (z0 - z)) * (
        np.exp(-((y - (ys + y0))**2 + x**2) / (2 * k**2 * (z0 - z)**2)) - 
        np.exp(-((y - yf)**2 + x**2) / (2 * k**2 * (z0 - z)**2))
    )
    
    return u_x, u_y, u_z

# ---- helpers (direct ports, minimal edits) ----
def _depth_parabolic(z_vec, Hw):
    z = np.asarray(z_vec).reshape(-1, 1)
    delta_w = (-6.0 / Hw) * z**2 + 6.0 * z
    return delta_w

def _depth_dmm(z_vec, Hw, C1, C2, C3):
    z = np.asarray(z_vec).reshape(-1, 1)
    delta_w = (
        C1 * (2*Hw - 2*z)
        + C2 * ((-6.0 / Hw) * z**2 + 6.0 * z)
        + C3 * (2.0 * z)
    )
    return delta_w

def _long_none(s_vec, Lwall, Hw, He_Hwratio):
    return np.ones((1, len(s_vec)), dtype=float)

def _long_mu(s_vec, Lwall, Hw, He_Hwratio):
    s = np.asarray(s_vec).reshape(1, -1)
    # width from Mu & Huang (kept as in your MATLAB)
    W = (Lwall/2.0) * (0.069*np.log((Hw/He_Hwratio)/Lwall) + 1.03)
    #W = max(W, 0.05*Lwall)
    return np.exp(-np.pi * (s / W) ** 2)


def _long_roboski(s_vec, Lwall, Hw, He_Hwratio):
    """
    Roboski & Finno (2006) erfc-type longitudinal reduction function.
    Vectorized equivalent of the MATLAB function:
        reduction = 1 - 0.5 .* erfc(num ./ den)
    """
    s = np.asarray(s_vec, dtype=float).ravel()  # ensure 1D
    a = 2.8
    c = 0.015 + 0.035 * np.log(He_Hwratio * Hw / Lwall)
    d = 0.15  + 0.035 * np.log(He_Hwratio * Hw / Lwall)
    num = a * ((Lwall / 2.0 - np.abs(s)) + Lwall * c)
    den = max((0.5 * Lwall - Lwall * d), 1e-9)  # avoid division by zero
    reduction = 1.0 - 0.5 * erfc(num / den)
    return reduction.reshape(1, -1)  # match MATLAB’s row-vector output

# Pinto/AF shaft solution: displacement from one equivalent spherical cavity
def _eq_shaft_3d_af(x, y, z, h, nu, a):
    R1 = np.sqrt(x**2 + y**2 + (z - h)**2)
    R2 = np.sqrt(x**2 + y**2 + (z + h)**2)

    fx = x * (1.0/R1**3 + (3 - 4*nu)/R2**3 - 6*z*(z + h)/R2**5)
    fy = y * (1.0/R1**3 + (3 - 4*nu)/R2**3 - 6*z*(z + h)/R2**5)
    fz = -1.0 * ((z - h)/R1**3 + 2*z/R2**3 - (3 - 4*nu)*(z + h)/R2**3 - 6*z*(z + h)**2/R2**5)

    # u = (1/3) * a^3 * f ; sign convention as in MATLAB (ux,uy negated)
    coeff = (a**3) / 3.0
    u_x = -coeff * fx
    u_y = -coeff * fy
    u_z =  coeff * fz
    return u_x, u_y, u_z

# Helper functions for DMM-based 3D greenfield model
def run_greenfield_3D_cloud_xyz(
    x, y, z, *,
    # geometry of the rectangular box (station half-lengths)
    Hw=19.0, L_x=9.5, L_y=32.0, He_Hwratio=1.0,
    # soil
    nu=0.499,
    # wall deflection shape (depth + longitudinal)
    switch_shape=5, C1=0.0, C2=1.0, C3=None,
    # equivalent-cavity amplitude per wall
    beta_CCS_wall_1=0.075/100, beta_CCS_wall_2=0.075/100,
    beta_CCS_wall_3=0.075/100, beta_CCS_wall_4=0.075/100,
    # discretization
    delta_z_cavities=1.0, delta_xyperimeter_cavities=2.5,
    # superposition type: 1=single wall mirrored, 2=four walls analytical, 3=four walls semi-analytical (with tapers)
    switch_solution_type=3,
    # zero-out inside box footprint? (greenfield outside only)
    zero_inside_box=True,
    debug=False,
):
    """
    DMM-based 3D greenfield at arbitrary (x,y,z) points.
    Inputs x, y, z must be same-shape arrays (1D ok). Returns u_x, u_y, u_z arrays.
    Requires helper funcs in this module: _depth_parabolic, _depth_dmm, _long_none,
    _long_mu, _long_roboski, _eq_shaft_3d_af.
    """
    if debug:
        print("\n===== DEBUG: run_greenfield_3D_cloud_xyz =====")
        print(f"Hw={Hw}, L_x={L_x}, L_y={L_y}, He_Hwratio={He_Hwratio}")
        print(f"nu={nu}")
        print(f"switch_shape={switch_shape}, switch_solution_type={switch_solution_type}")
        print(f"C1={C1}, C2={C2}, C3={C3}")
        print("betas:",
              beta_CCS_wall_1, beta_CCS_wall_2,
              beta_CCS_wall_3, beta_CCS_wall_4)
        print(f"delta_z_cavities={delta_z_cavities}, "
              f"delta_xyperimeter_cavities={delta_xyperimeter_cavities}")
   
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    if not (x.shape == y.shape == z.shape):
        raise ValueError("x, y, z must have the same shape.")
    Np = x.size

    # depth weights
    if C3 is None:
        C3 = 1.0 - C1 - C2
    if switch_shape in (5, 50, 51) and abs((C1+C2+C3) - 1.0) > 1e-6:
        raise ValueError("For DMM depth (5/50/51), C1+C2+C3 must equal 1.")

    # choose depth function
    if switch_shape in (3, 30, 31):
        depth_fun = lambda zz: _depth_parabolic(zz, Hw)
    elif switch_shape in (5, 50, 51):
        depth_fun = lambda zz: _depth_dmm(zz, Hw, C1, C2, C3)
    else:
        raise ValueError(f"Unknown switch_shape={switch_shape}")

    # choose longitudinal reduction
    if switch_shape in (3, 5):
        long_fun = _long_none
    elif switch_shape in (30, 50):
        long_fun = _long_mu
    elif switch_shape in (31, 51):
        long_fun = _long_roboski
    else:
        raise ValueError(f"Unknown switch_shape={switch_shape}")

    # perimeter discretization (four walls)
    n1 = max(1, int(round(2*L_x / delta_xyperimeter_cavities)))
    n2 = max(1, int(round(2*L_y / delta_xyperimeter_cavities)))
    dx1 = (2*L_x) / n1
    dx2 = (2*L_y) / n2

    # wall midpoints
    xmid_w1 = np.linspace(-L_x + dx1/2,  L_x - dx1/2, n1); ymid_w1 = -L_y*np.ones(n1)
    ymid_w2 = np.linspace(-L_y + dx2/2,  L_y - dx2/2, n2); xmid_w2 =  L_x*np.ones(n2)
    xmid_w3 = np.linspace( L_x - dx1/2, -L_x + dx1/2, n1); ymid_w3 =  L_y*np.ones(n1)
    ymid_w4 = np.linspace( L_y - dx2/2, -L_y + dx2/2, n2); xmid_w4 = -L_x*np.ones(n2)

    Xc_vec = np.concatenate([xmid_w1, xmid_w2, xmid_w3, xmid_w4])
    Yc_vec = np.concatenate([ymid_w1, ymid_w2, ymid_w3, ymid_w4])

    idx_w1 = np.arange(0, n1)
    idx_w2 = np.arange(n1, n1+n2)
    idx_w3 = np.arange(n1+n2, 2*n1+n2)
    idx_w4 = np.arange(2*n1+n2, 2*n1+2*n2)

    Lwall_1 = 2*L_x;  Lwall_2 = 2*L_y;  Lwall_3 = 2*L_x;  Lwall_4 = 2*L_y

    # vertical discretization (cavity mid-depths)
    Nlayers = max(1, int(round(Hw / delta_z_cavities)))
    delta_z_cavities = Hw / Nlayers  # adjust to fit exactly
    z_wall_discr = np.linspace(0.0, Hw, Nlayers + 1)
    z_cav = 0.5 * (z_wall_discr[:-1] + z_wall_discr[1:])
    Nd = z_cav.size

    if debug:
        print(f"Nlayers={Nlayers}, final delta_z_cavities={delta_z_cavities}")
        print("z_wall_discr (first few):", z_wall_discr[:5])
        print("z_cav (first few):", z_cav[:5])

    # amplitude matrix: depth × longitudinal along each wall
    delta_cav = depth_fun(z_cav.reshape(-1, 1))  # Nd×1
    R1 = long_fun(xmid_w1, Lwall_1, Hw, He_Hwratio)  # 1×n1
    R2 = long_fun(ymid_w2, Lwall_2, Hw, He_Hwratio)  # 1×n2
    R3 = long_fun(xmid_w3, Lwall_3, Hw, He_Hwratio)  # 1×n1
    R4 = long_fun(ymid_w4, Lwall_4, Hw, He_Hwratio)  # 1×n2

    A1 = beta_CCS_wall_1 * (delta_cav @ R1)  # Nd×n1
    A2 = beta_CCS_wall_2 * (delta_cav @ R2)  # Nd×n2
    A3 = beta_CCS_wall_3 * (delta_cav @ R3)  # Nd×n1
    A4 = beta_CCS_wall_4 * (delta_cav @ R4)  # Nd×n2
    A = np.concatenate([A1, A2, A3, A4], axis=1)     # Nd×Ns
    
    if debug:
        # Integrate along depth to approximate "volume loss per run meter"
        VLW_w1_py = np.trapezoid(A1, z_cav, axis=0)   # length n1
        VLW_w2_py = np.trapezoid(A2, z_cav, axis=0)   # length n2
        VLW_w3_py = np.trapezoid(A3, z_cav, axis=0)
        VLW_w4_py = np.trapezoid(A4, z_cav, axis=0)

        beta_run_w1_py = VLW_w1_py / (Hw**2)
        beta_run_w2_py = VLW_w2_py / (Hw**2)
        beta_run_w3_py = VLW_w3_py / (Hw**2)
        beta_run_w4_py = VLW_w4_py / (Hw**2)

        print("beta_runmeter_wall_2 (Python, first 5):",
              beta_run_w2_py[:5])

    # symmetry factor & segment lengths
    if switch_solution_type == 1:
        m_sym = 2.0
    elif switch_solution_type == 2:
        m_sym = 1.0
    elif switch_solution_type == 3:
        m_sym = 2.0
    else:
        raise ValueError("switch_solution_type must be 1, 2, or 3")

    seg_len = np.concatenate([
        np.full(n1, dx1),
        np.full(n2, dx2),
        np.full(n1, dx1),
        np.full(n2, dx2),
    ])[None, :]  # 1×Ns

    Nd, Ns = A.shape
    vol_elem = m_sym * A * (delta_z_cavities * seg_len)    # Nd×Ns
    X = (3.0 / (4.0*np.pi)) * vol_elem # Allow for negative volumes to be passed through
    a_cavity = np.sign(X) * np.abs(X) ** (1.0 / 3.0) 

    #vol_elem = np.maximum(vol_elem, 0.0)
    #a_cavity = ((3.0 / (4.0*np.pi)) * vol_elem) ** (1.0/3.0)  # Nd×Ns

    if debug:
        print("vol_elem min/max:", float(vol_elem.min()), float(vol_elem.max()))
        print("a_cavity min/max:", float(a_cavity.min()), float(a_cavity.max()))

    # per-point tapers for semi-analytical option (type 3)
    if switch_solution_type == 3:
        F1  = np.ones(Np, dtype=float)
        F2  = np.ones(Np, dtype=float)
        F3  = np.ones(Np, dtype=float)
        F4  = np.ones(Np, dtype=float)

        # W1 (bottom, Y=-Ly): 1 at Y=-Ly → 0 at Y=+Ly, only for y >= -Ly
        mask_w1 = (y >= -L_y)
        F1[mask_w1] = np.maximum(0.0, 0.5 - 0.5*(y[mask_w1] / L_y))  # 1 at y=-Ly → 0 at y=+Ly
        # W2 (right, X=+Lx): 1 at X=+Lx → 0 at X=-Lx, only for x <= +Lx
        mask_w2 = (x <= L_x)
        F2[mask_w2] = np.maximum(0.0, 0.5 + 0.5*(x[mask_w2] / L_x))  # 1 at x=+Lx → 0 at x=-Lx
        # W3 (top, Y=+Ly): 1 at Y=+Ly → 0 at Y=-Ly, only for y <= +Ly
        mask_w3 = (y <= L_y)
        F3[mask_w3] = np.maximum(0.0, 0.5 + 0.5*(y[mask_w3] / L_y))  # 1 at y=+Ly → 0 at y=-Ly
        # W4 (left, X=-Lx): 1 at X=-Lx
        mask_w4 = (x >= -L_x)
        F4[mask_w4] = np.maximum(0.0, 0.5 - 0.5*(x[mask_w4] / L_x))  # 1 at x=-Lx → 0 at x=+Lx
    else:
        F1 = F2 = F3 = F4 = np.ones(Np)

    # accumulate contributions
    ux_all = np.zeros(Np, dtype=float)
    uy_all = np.zeros(Np, dtype=float)
    uz_all = np.zeros(Np, dtype=float)

    index_sets = [idx_w1, idx_w2, idx_w3, idx_w4]
    F_map = {1: F1, 2: F2, 3: F3, 4: F4}

    for w, idxs in enumerate(index_sets, start=1):
        Fw = F_map[w]
        for j in idxs:
            dxj = Xc_vec[j]; dyj = Yc_vec[j]
            xloc = x - dxj
            yloc = y - dyj
            zloc = z
            for i in range(Nd):
                a = a_cavity[i, j]
                #if a <= 0.0: Allow for negative cavities
                 #   continue
                h = z_cav[i]
                ux, uy, uz = _eq_shaft_3d_af(xloc, yloc, zloc, h, nu, a)
                ux_all += Fw * ux
                uy_all += Fw * uy
                uz_all += Fw * uz

    if zero_inside_box:
        inside = (x >= -L_x) & (x <= L_x) & (np.abs(y) <= L_y)
        ux_all[inside] = 0.0
        uy_all[inside] = 0.0
        uz_all[inside] = 0.0

    return (np.ascontiguousarray(ux_all, dtype=np.float64),
            np.ascontiguousarray(uy_all, dtype=np.float64),
            np.ascontiguousarray(uz_all, dtype=np.float64))


# Approach 2 for C1 C2 and C3 variation along the walls

def run_greenfield_3D_cloud_xyz_approach2(
        x,y,z, *,
        # geometry of the rectangular box (station half-lengths)
        Hw=19.0, L_x=9.5, L_y=32.0, He_Hwratio=1.0,
        # soil
        nu=0.499,
        # wall deflection shape (DMM + Mu/Roboski)
        switch_shape = 50,
        # equivalent-cavity amplitude per wall (beta at centre line)
        beta_CCS_wall_1=0.075/100, beta_CCS_wall_2=0.075/100,
        beta_CCS_wall_3=0.075/100, beta_CCS_wall_4=0.075/100,

        #Discretization
        delta_z_cavities=1.0, delta_xyperimeter_cavities=2.5,
        # superposition type: 1=single wall mirrored, 2=four walls analytical, 3=four walls semi-analytical (with tapers)
        switch_solution_type=3,
        # zero-out inside box footprint? (greenfield outside only)
        zero_inside_box=True,
        # Anchor Inputs for Approach 2
        y_anchor_long = None, 
        x_anchor_short = None, 
        C1_anchor_w1 = None, C2_anchor_w1 = None, C3_anchor_w1 = None,
        C1_anchor_w2 = None, C2_anchor_w2 = None, C3_anchor_w2 = None,
        C1_anchor_w3 = None, C2_anchor_w3 = None, C3_anchor_w3 = None,
        C1_anchor_w4 = None, C2_anchor_w4 = None, C3_anchor_w4 = None,
):
    
    """
    DMM-based 3D greenfield at arbitrary (x,y,z) points.
    C1, C2 C2 variation along the walls based on user defined anchor points.
    -Same geometry as run_greenfield_3D_cloud_xyz
    - switch_shape must be:
        50 : DMM + Mu
        51 : DMM + Roboski
    - C1, C2 and C3 along each wall are defined by anchor points and linearly interpolated to the cavity midpoints on each wall. 
    Returns u_x, u_y, u_z arrays.
    """

    # -- Basic checks and input shapes ---
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    if not (x.shape == y.shape == z.shape):
        raise ValueError("x, y, z must have the same shape.")
    Np = x.size

    # -- Restrict to DMM + Mu/Roboski ---
    if switch_shape == 50:
        long_fun = _long_mu
    elif switch_shape == 51:
        long_fun = _long_roboski
    else:
        raise ValueError("For Approach 2 (run_greenfield_3D_cloud_xyz_approach2) currently supports switch_shape= 50 (Mu) or 51 (Roboski).")
    
    Lwall_1 = 2*L_x
    Lwall_2 = 2*L_y
    Lwall_3 = 2*L_x
    Lwall_4 = 2*L_y

    # -- Default for anchor positions ---

    if y_anchor_long is None:
        y_anchor_long = np.array(
            [0.0, 0.25*Lwall_2, 0.5*Lwall_2, 0.75*Lwall_2, Lwall_2]
            )
    else:
        y_anchor_long = np.asarray(y_anchor_long, dtype=float)

    if x_anchor_short is None:
        x_anchor_short = np.array(
            [0.0,  0.5*Lwall_1,  Lwall_1]
            )
    else:
        x_anchor_short = np.asarray(x_anchor_short, dtype=float)

    # helper for default anchors
    def _default(arr,fallback):
        return np.asarray(arr, dtype=float) if arr is not None else np.asarray(fallback, dtype=float)
    
    # -- Default Values for C1, C2, C3 anchors (correspond to KK) ---
    C1_anchor_w2 = _default(C1_anchor_w2, [0.33, -0.5, 0.3439, 0.26, 0.33])
    C2_anchor_w2 = _default(C2_anchor_w2, [0.33,  0.26, 0.1769, -0.09, 0.33])
    C3_anchor_w2 = _default(C3_anchor_w2, [0.34,  1.24, 0.4792,  0.83, 0.34])

    C1_anchor_w4 = _default(C1_anchor_w4, [0.33, -0.5, 0.3439, 0.26, 0.33])
    C2_anchor_w4 = _default(C2_anchor_w4, [0.33,  0.26, 0.1769, -0.09, 0.33])
    C3_anchor_w4 = _default(C3_anchor_w4, [0.34,  1.24, 0.4792,  0.83, 0.34])

    C1_anchor_w1 = _default(C1_anchor_w1, [0.33, 0.40, 0.33])
    C2_anchor_w1 = _default(C2_anchor_w1, [0.33, 0.20, 0.33])
    C3_anchor_w1 = _default(C3_anchor_w1, [0.34, 0.40, 0.34])

    C1_anchor_w3 = _default(C1_anchor_w3, [0.33, 0.40, 0.33])
    C2_anchor_w3 = _default(C2_anchor_w3, [0.33, 0.20, 0.33])
    C3_anchor_w3 = _default(C3_anchor_w3, [0.34, 0.40, 0.34])

    # -- Vertical Discretization in z ---
    z_wall_discr = np.arange(0.0, Hw + delta_z_cavities, delta_z_cavities)
    if z_wall_discr.size < 2:
        z_wall_discr = np.array([0.0, Hw], dtype=float)
    z_cav = 0.5 * (z_wall_discr[:-1] + z_wall_discr[1:])
    Nd_c = z_cav.size

    # -- Perimeter discretization (four walls) ---
    n1 = max(1, int(round(2*L_x / delta_xyperimeter_cavities)))
    n2 = max(1, int(round(2*L_y / delta_xyperimeter_cavities)))
    dx1 = (2*L_x) / n1
    dx2 = (2*L_y) / n2
    # wall midpoints
    xmid_w1 = np.linspace(-L_x + dx1/2,  L_x - dx1/2, n1); ymid_w1 = -L_y*np.ones(n1)
    ymid_w2 = np.linspace(-L_y + dx2/2,  L_y - dx2/2, n2); xmid_w2 =  L_x*np.ones(n2)
    xmid_w3 = np.linspace( L_x - dx1/2, -L_x + dx1/2, n1); ymid_w3 =  L_y*np.ones(n1)
    ymid_w4 = np.linspace( L_y - dx2/2, -L_y + dx2/2, n2); xmid_w4 = -L_x*np.ones(n2)

    Xc_vec = np.concatenate([xmid_w1, xmid_w2, xmid_w3, xmid_w4])
    Yc_vec = np.concatenate([ymid_w1, ymid_w2, ymid_w3, ymid_w4])
    idx_w1 = np.arange(0, n1)
    idx_w2 = np.arange(n1, n1+n2)
    idx_w3 = np.arange(n1+n2, 2*n1+n2)
    idx_w4 = np.arange(2*n1+n2, 2*n1+2*n2)

    Ns = Xc_vec.size # total number of perimeter segments/stacks

    # -- Along wall coordinates --

    s_w1 = xmid_w1  # [-Lx+dx/2, ..., Lx-dx/2]
    s_w2 = ymid_w2  # [-Ly+dy/2, ..., Ly-dy/2]
    s_w3 = xmid_w3  # [Lx-dx/2, ..., -Lx+dx/2]
    s_w4 = ymid_w4  # [Ly-dy/2, ..., -Ly+dy/2]

    # -- Longitudinal reduction R1, R2, R3 and R4 --
    R1 = long_fun(s_w1, Lwall_1, Hw, He_Hwratio)  # 1×n1
    R2 = long_fun(s_w2, Lwall_2, Hw, He_Hwratio)  # 1×n2
    R3 = long_fun(s_w3, Lwall_3, Hw, He_Hwratio)  # 1×n1
    R4 = long_fun(s_w4, Lwall_4, Hw, He_Hwratio)  # 1×n2

    # reshape to 1D arrays for interpolation

    R1 = R1.reshape(-1)
    R2 = R2.reshape(-1)
    R3 = R3.reshape(-1)
    R4 = R4.reshape(-1)

    # -- Interpolate C1, C2, C3 along each wall to cavity midpoints --
    # Long walls: s + L_y → [0, 2L_y]; Short walls: s + L_x → [0, 2L_x]
    C1_prof_w2 = np.interp(s_w2 + L_y, y_anchor_long, C1_anchor_w2)
    C2_prof_w2 = np.interp(s_w2 + L_y, y_anchor_long, C2_anchor_w2)
    C3_prof_w2 = np.interp(s_w2 + L_y, y_anchor_long, C3_anchor_w2)

    C1_prof_w4 = np.interp(s_w4 + L_y, y_anchor_long, C1_anchor_w4)
    C2_prof_w4 = np.interp(s_w4 + L_y, y_anchor_long, C2_anchor_w4)
    C3_prof_w4 = np.interp(s_w4 + L_y, y_anchor_long, C3_anchor_w4)

    C1_prof_w1 = np.interp(s_w1 + L_x, x_anchor_short, C1_anchor_w1)
    C2_prof_w1 = np.interp(s_w1 + L_x, x_anchor_short, C2_anchor_w1)
    C3_prof_w1 = np.interp(s_w1 + L_x, x_anchor_short, C3_anchor_w1)

    C1_prof_w3 = np.interp(s_w3 + L_x, x_anchor_short, C1_anchor_w3)
    C2_prof_w3 = np.interp(s_w3 + L_x, x_anchor_short, C2_anchor_w3)
    C3_prof_w3 = np.interp(s_w3 + L_x, x_anchor_short, C3_anchor_w3)

    # -- Build delta_cavity for all walls --

    Delta_cav_all = np.zeros((Nd_c, Ns), dtype=float)

    # Wall 1

    for j, j_global in enumerate(idx_w1):
        delta_cav_j = _depth_dmm(z_cav, Hw, C1_prof_w1[j], C2_prof_w1[j], C3_prof_w1[j]).reshape(-1)
        Delta_cav_all[:, j_global] = (beta_CCS_wall_1 * R1[j]) * delta_cav_j

    # Wall 2

    for j, j_global in enumerate(idx_w2):
        jj = j  # local index for wall 2
        delta_cav_j = _depth_dmm(z_cav, Hw, C1_prof_w2[jj], C2_prof_w2[jj], C3_prof_w2[jj]).reshape(-1)
        Delta_cav_all[:, j_global] = (beta_CCS_wall_2 * R2[jj]) * delta_cav_j

    # Wall 3
    for j, j_global in enumerate(idx_w3):
        jj = j  # local index for wall 3
        delta_cav_j = _depth_dmm(z_cav, Hw, C1_prof_w3[jj], C2_prof_w3[jj], C3_prof_w3[jj]).reshape(-1)
        Delta_cav_all[:, j_global] = (beta_CCS_wall_3 * R3[jj]) * delta_cav_j

    # Wall 4
    for j, j_global in enumerate(idx_w4):
        jj = j  # local index for wall 4
        delta_cav_j = _depth_dmm(z_cav, Hw, C1_prof_w4[jj], C2_prof_w4[jj], C3_prof_w4[jj]).reshape(-1)
        Delta_cav_all[:, j_global] = (beta_CCS_wall_4 * R4[jj]) * delta_cav_j
     
    # -- Symmetry factor & segment lengths --
    if switch_solution_type == 1:
        m_sym = 2.0
    elif switch_solution_type == 2:
        m_sym = 1.0
    elif switch_solution_type == 3:
        m_sym = 2.0
    else:
        raise ValueError("switch_solution_type must be 1, 2, or 3")
    
    seg_len = np.concatenate([
        np.full(n1, dx1),
        np.full(n2, dx2),
        np.full(n1, dx1),
        np.full(n2, dx2),
    ])

    # -- Volume element and equivalent cavity radius (Allow for negative cavities) --

    vol_elem = m_sym * Delta_cav_all * (delta_z_cavities * seg_len)     # Nd×Ns
    Xvol = (3.0 / (4.0*np.pi)) * vol_elem
    a_cavity = np.sign(Xvol) * np.abs(Xvol) ** (1.0 / 3.0)  # Nd×Ns

    # -- Per-point tapers for semi-analytical option (type 3) --
    if switch_solution_type == 3:
        F1  = np.ones(Np, dtype=float)
        F2  = np.ones(Np, dtype=float)
        F3  = np.ones(Np, dtype=float)
        F4  = np.ones(Np, dtype=float)

        # W1 (bottom, Y=-Ly): 1 at Y=-Ly → 0 at Y=+Ly, only for y >= -Ly
        mask_w1 = (y >= -L_y)
        F1[mask_w1] = np.maximum(0.0, 0.5 - 0.5*(y[mask_w1] / L_y))  # 1 at y=-Ly → 0 at y=+Ly
        # W2 (right, X=+Lx): 1 at X=+Lx → 0 at X=-Lx, only for x <= +Lx
        mask_w2 = (x <= L_x)
        F2[mask_w2] = np.maximum(0.0, 0.5 + 0.5*(x[mask_w2] / L_x))  # 1 at x=+Lx → 0 at x=-Lx
        # W3 (top, Y=+Ly): 1 at Y=+Ly → 0 at Y=-Ly, only for y <= +Ly
        mask_w3 = (y <= L_y)
        F3[mask_w3] = np.maximum(0.0, 0.5 + 0.5*(y[mask_w3] / L_y))  # 1 at y=+Ly → 0 at y=-Ly
        # W4 (left, X=-Lx): 1 at X=-Lx
        mask_w4 = (x >= -L_x)
        F4[mask_w4] = np.maximum(0.0, 0.5 - 0.5*(x[mask_w4] / L_x))  # 1 at x=-Lx → 0 at x=+Lx
        
    else:
        F1 = F2 = F3 = F4 = np.ones(Np)

    # -- Accumulate contributions --
    ux_all = np.zeros(Np, dtype=float)
    uy_all = np.zeros(Np, dtype=float)
    uz_all = np.zeros(Np, dtype=float)

    index_sets = [idx_w1, idx_w2, idx_w3, idx_w4]
    F_map = {1: F1, 2: F2, 3: F3, 4: F4}

    for w, idxs in enumerate(index_sets, start=1):
        Fw = F_map[w]
        for j in idxs:
            dxj = Xc_vec[j]; dyj = Yc_vec[j]
            xloc = x - dxj
            yloc = y - dyj
            zloc = z
            for i in range(Nd_c):
                a = a_cavity[i, j]
                #if a <= 0.0: Allow for negative cavities
                 #   continue
                h = z_cav[i]
                ux, uy, uz = _eq_shaft_3d_af(xloc, yloc, zloc, h, nu, a)
                ux_all += Fw * ux
                uy_all += Fw * uy
                uz_all += Fw * uz

    if zero_inside_box:
        inside = (x >= -L_x) & (x <= L_x) & (np.abs(y) <= L_y)
        ux_all[inside] = 0.0
        uy_all[inside] = 0.0
        uz_all[inside] = 0.0

    return (np.ascontiguousarray(ux_all, dtype=np.float64),
            np.ascontiguousarray(uy_all, dtype=np.float64),
            np.ascontiguousarray(uz_all, dtype=np.float64))