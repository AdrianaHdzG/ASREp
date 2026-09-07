# %% Imports
import sys
import os
import numpy as np

# Path setup: find ASREp package and Input.py regardless of working directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))               # Input.py
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..'))                                   # ASREp package
sys.path.insert(0, r"C:\ASREp")
from Input import *
from ASREp.ground_deformation_models import wall_deflection
from ASREp.ASRE_Timoshenko_model import ASRE_Timoshenko_model

# %% GREENFIELD — FULL 2D GRID (for contour plots)
print('--- Computing greenfield on full 2D grid ---')

(uh_grid, uz_grid,
 wall_deflection_profile, VLW, max_wall_defl,
 check_ratio, VLS, ratio_vol) = wall_deflection(
    retaining_wall_depth, avg_wall_displacement, wall_deflection_shape,
    soil_ovalization, volumetric,
    x_coords=x_grid, z_coords=z_grid,
    C1_val=C1_val, C2_val=C2_val
)

# Depth vector matching wall_deflection_profile (internal step = Hwall/100)
deltaH    = retaining_wall_depth / 100
depth_vec = np.arange(0, retaining_wall_depth + deltaH, deltaH)

print(f'  Max wall deflection      : {max_wall_defl * 1000:.2f} mm')
print(f'  Wall volume loss VLW     : {VLW:.4f} m²')
print(f'  Surface/wall vol ratio   : {ratio_vol:.3f}')
print(f'  Integration check ratio  : {check_ratio:.4f}  (should be ~1)')

# %% GREENFIELD — BUILDING LINE (at foundation depth)
print('\n--- Computing greenfield at building line ---')

(uh_building_2d, uz_building_2d, _, _, _, _, _, _) = wall_deflection(
    retaining_wall_depth, avg_wall_displacement, wall_deflection_shape,
    soil_ovalization, volumetric,
    x_coords=x_coords_building, z_coords=z_coords_building,
    C1_val=C1_val, C2_val=C2_val
)

# Flatten from (1, num_nodes) to (num_nodes,)
Sx = uh_building_2d.flatten().astype(np.float64)   # [m] horizontal, positive away from wall
Sz = uz_building_2d.flatten().astype(np.float64)   # [m] vertical,   positive downward (geotechnical)

# %% GREENFIELD STRAIN ANALYSIS (Boscardin & Cording 1989)
print('\n--- Greenfield strain analysis ---')

length_element    = building_length / num_elements
tilt              = (Sz[-1] - Sz[0]) / building_length    # [rad] rigid body tilt
x_elem_centres    = 0.5 * (x_coords_building[:-1] + x_coords_building[1:])

S       = np.diff(Sz) / length_element                    # slope of settlement profile [rad]
beta_d  = S - tilt                                        # angular distortion, tilt removed [rad]
eps_h   = np.diff(Sx) / length_element                    # horizontal strain [-]

# Maximum principal tensile strain from eps_h and beta_d (Mohr's circle, eps_v = 0)
eps_t_gf = 0.5 * eps_h + np.sqrt((0.5 * eps_h) ** 2 + (0.5 * beta_d) ** 2)
eps_t_gf = np.maximum(eps_t_gf, 0.0)

max_eps_t_gf = float(np.max(eps_t_gf))
print(f'  Max angular distortion   : {np.max(np.abs(beta_d)) * 1000:.4f} ‰')
print(f'  Max horizontal strain    : {np.max(np.abs(eps_h)) * 100:.4f} %')
print(f'  Max tensile strain (GF)  : {max_eps_t_gf * 100:.4f} %')

# %% TIMOSHENKO SSI MODEL
print('\n--- Running Timoshenko SSI model ---')

# Sign convention (matches MastersThesis main.py):
# wall_deflection returns Sz positive-downward (geotechnical convention)
# Timoshenko expects upward-positive (FEM convention)  →  negate Sz
# UPWARDS POSITIVE, TOWARDS RIGHT POSITIVE AND COUNTERCLOCKWISE POSITIVE
dispL = Sx                                      # [m] longitudinal (x), positive away from wall
dispT = np.zeros(num_nodes, dtype=np.float64)   # [m] transverse  (y), zero in 2D
dispV = -Sz                                     # [m] vertical, FEM convention  # STD: FEM CONVENTION

beamY = np.zeros_like(x_coords_building, dtype=np.float64)
beamZ = np.zeros_like(x_coords_building, dtype=np.float64)

model = ASRE_Timoshenko_model(
    num_nodes, x_coords_building, beamY, beamZ,
    building_height, building_width,
    d_NA=d_NA, solver='elastic'
)
model.set_beam_properties(Eb, EoverG, q_foot, d_NA=d_NA)
model.set_soil_properties(Es_isotropic, soil_poisson, mu_int)

model.run_model(dispL, dispT, dispV, output='strain+disp+force')

damage_category, max_eps_t_ssi = model.categorize_damage()

print(f'  Max tensile strain (SSI) : {max_eps_t_ssi:.4f} %')
print(f'  Damage category (SSI)    : {damage_category}  '
      f'(0=negligible -> 5=very severe)')

# %% SUMMARY
print('\n========== SUMMARY ==========')
print(f'  Wall deflection shape    : {wall_deflection_shape}  '
      f'(0=Uniform,1=Cantilever,2=Parabola,3=Composite,4=Kick-in,5=Custom)')
print(f'  Avg wall displacement    : {avg_wall_displacement * 100:.3f} % of Hwall')
print(f'  Max wall deflection      : {max_wall_defl * 1000:.2f} mm')
print(f'  GF  max tensile strain   : {max_eps_t_gf * 100:.4f} %')
print(f'  SSI max tensile strain   : {max_eps_t_ssi:.4f} %')
print(f'  Damage category (SSI)    : {damage_category}')


# --------------------------------------------------
# quoFEM output
# --------------------------------------------------

with open("results.out", "w") as f:
    f.write(
        "{} {} {}".format(
            max_eps_t_gf,
            max_eps_t_ssi,
            damage_category
        )
    )

print("results.out written successfully")