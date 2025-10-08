%% run_greenfield_3D_DMM.m
% -------------------------------------------------------------------------
% Purpose:
%   Compute 3D greenfield displacements induced by excavation, based on
%   user-defined retaining wall deflection shapes.
%
% Key Features:
%   1. Supports multiple wall deflection modes:
%        - 3  : Parabolic shape (with/without degradation)
%        - 30 : Parabolic with Gaussian longitudinal reduction (Mu & Huang 2016)
%        - 31 : Parabolic with erfc-type longitudinal reduction (Roboski & Finno 2006)
%        - 5  : Custom linear combination of cantilever, parabolic, and kick-in
%   2. Allows specification of wall displacement level  beta and
%      variation along walls (longitudinal direction).
%   3. Discretizes excavation cavity into vertical slices (depth) and
%      perimeter segments (plan view).
%   4. Computes cavity area contributions, wall displacement fields,
%      and beta/ beta_CSS ratios for verification.
%
% Inputs (set in script):
%   - switch_shape        : integer selecting wall deflection shape
%   - beta_CCS_wall_i     : control  beta for each wall (1–4)
%   - Hw, Hr, Lx, Ly      : geometric parameters (wall height, reference depth,
%                            half-lengths of excavation box)
%   - delta_z_cavities,
%     delta_xyperimeter   : discretization step sizes (depth & perimeter)
%
% Outputs (variables/plots):
%   - a_cavity_matrix_v*  : cavity displacement contributions per slice
%   - VLW_runmeter_wall_* : integrated wall deflections per unit length
%   - beta/ beta_CSS diagnostic plots (figure 97, 98)
%
% Notes:
%   - Discretization along walls is midpointed to ensure symmetry.
%   - Longitudinal reductions (30/31) are centered around s = 0 for each wall.
%   - Use figure(97)/(98) to check symmetry of β/β_CSS profiles.
%
% Author: Adriana Hernandez
% Created: [Date]
% Based on: initial 3DShaft framework (2025), adapted for 3D wall
% deflection created by Andrea Franza
% _v1: Station_calculation_AFv4_Eq_shaft_3d_AF_v2
% _v2: Updates discretization and perimeter location of cavities. 
% -------------------------------------------------------------------------

%close all;
%clear all;
function out = run_greenfield_3D_DMM(params)
    % Run the 3D greenfield DMM with optional overrides via params struct.
    % Example:
    %   p = struct('Hw',19,'L_x',9.5,'L_y',32,'switch_shape',5);
    %   out = run_greenfield_3D_DMM(p);
                          % ← ensures out is always defined
    if nargin < 1, params = struct(); end
    out = struct(); 
    close all; clc;  % keep it light; avoid clearvars/clear all inside functions
    
    %% ---------------- Inputs (guarded defaults) ----------------
    nu            = getf(params,'nu',            0.499);   % Poisson's ratio
    Hw            = getf(params,'Hw',            19);      % wall depth [m]
    He_Hwratio    = getf(params,'He_Hwratio',    1);       % He/Hw ratio if used
    L_x           = getf(params,'L_x',           9.5);     % half-width of box [m]
    L_y           = getf(params,'L_y',           32);      % half-height of box [m]
    
    beta_CCS_wall_1 = getf(params,'beta_CCS_wall_1', 0.075/100);
    beta_CCS_wall_2 = getf(params,'beta_CCS_wall_2', 0.075/100);
    beta_CCS_wall_3 = getf(params,'beta_CCS_wall_3', 0.075/100);
    beta_CCS_wall_4 = getf(params,'beta_CCS_wall_4', 0.075/100);
    
    
    num_nodes       = getf(params,'num_nodes',   101);     % number of evaluation points along the line
    building_offset = getf(params,'building_offset',   11);    % [m] distance from the wall face
    length_beam     = getf(params,'length_beam',   12);    % [m] length along x
    y0              = getf(params,'y0',   0);     % [m] fixed y for the line
    z0              =getf(params,'z0',   0);      % [m] fixed z for the line (0=surface; e.g. 2.2 for foundation)     
    
    
    switch_shape   = getf(params,'switch_shape',   5);     % 3,30,31,5,50,51 (your modes)
    C1             = getf(params,'C1',             0);     % only used if switch_shape=5/50/51
    C2             = getf(params,'C2',             1);
    C3             = getf(params,'C3',             1-C1-C2);

    if ismember(switch_shape, [5 50 51])
        tol  = 1e-6;
        if abs((C1+C2+C3) - 1) > tol
            error('C1+C2+C3 must equal 1 (got %.6g). Fix your inputs.', C1+C2+C3);
        end
    end
    
    % Wall discretization
    delta_z_cavities = getf(params,'delta_z_cavities', Hw/19);  % depth step [m]
    % Perimeter discretization
    delta_xyperimeter_cavities = getf(params,'delta_xyperimeter_cavities', 2.5);
    
    % Field grid spacing (used later when building Xsection/Ysection/Zsection)
    spacing_x = getf(params,'spacing_x', Hw*0.05);
    spacing_y = getf(params,'spacing_y', Hw*0.05);
    spacing_z = getf(params,'spacing_z', Hw*0.05);
    
    % Output section choice and plotting/solution switches
    switch_outputlocation   = getf(params,'switch_outputlocation',   3); % 1=xz, 3=xy, 4 =buildingcoords
    switch_solution_type    = getf(params,'switch_solution_type',    3); % 1,2,3 per your script
    
    % Make C’s visible if shape=5/50/51
    if ismember(switch_shape,[5 50 51]) 
        % expose C1..C3 to original code that references them
    end

    make_plots_checks = getf(params,'make_plots_checks', false);
    make_plots_results = getf(params,'make_plots_results', false);


    %% Building up coordinates at which movement computed
if switch_outputlocation==1
    %%xz plane (vertical)
    x_vect = L_x+spacing_x : spacing_x : 6*Hw;
    y_vect = 0;
    z_vect = 0 : spacing_z : 2*Hw;

elseif switch_outputlocation==3
    % xy plane (horizontal)
    x_vect = -5*Hw : spacing_x : 5*Hw;
    y_vect = -5*Hw : spacing_y : 5*Hw;
    z_vect = 0;

elseif switch_outputlocation==4
    % ---- 1D line along x, at fixed y=z ----
    % Global x along the building line (outside the RIGHT wall at x=+L_x):
    %   absolute x = (wall face) + (offset) + [0..length_beam]
    x_vect = (L_x + building_offset) + linspace(0, length_beam, num_nodes);
    y_vect = y0;     % scalar
    z_vect = z0;     % scalar

elseif switch_outputlocation==100
    % --- ASRE3D building nodes (direct) ---
    % Expect the caller to pass coords via params.x_nodes / y_nodes / z_nodes
    assert(isfield(params,'x_nodes') && isfield(params,'y_nodes') && isfield(params,'z_nodes'), ...
           'Mode 100 needs params.x_nodes, params.y_nodes, params.z_nodes');
    x_vect = params.x_nodes(:);
    y_vect = params.y_nodes(:);
    z_vect = params.z_nodes(:);

else
    error('switch_outputlocation=%d not implemented.', switch_outputlocation);
end


if switch_outputlocation==1
    %%xz plane (vertical)     
    [Xsection,Zsection] = meshgrid(x_vect,z_vect);
    Ysection=0;
   
elseif  switch_outputlocation==3
    %xy plane (horizontal)
    [Xsection,Ysection] = meshgrid(x_vect,y_vect);
    Zsection=0;

elseif switch_outputlocation==100
    % --- ASRE3D node cloud: 1-to-1 points, no meshgrid ---
    % x_vect, y_vect, z_vect must be same-length vectors of node coordinates.
    if ~(numel(x_vect)==numel(y_vect) && numel(y_vect)==numel(z_vect))
        error('Mode 100: x_nodes, y_nodes, z_nodes must have the same length.');
    end
    Xsection = x_vect(:);
    Ysection = y_vect(:);
    Zsection = z_vect(:);


else
    [Xsection,Ysection,Zsection] = meshgrid(x_vect,y_vect,z_vect);
    Xsection=Xsection(:);
    Ysection=Ysection(:);
    Zsection=Zsection(:);
end


%% Basic wall mode
    beta_CCS_wall_all = [beta_CCS_wall_1 , beta_CCS_wall_2 , beta_CCS_wall_3 , beta_CCS_wall_4];
    beta_CCS_avg = mean(beta_CCS_wall_all);
    
    z_wall_discr_vec = 0:delta_z_cavities:Hw; %depth at which wall volume loss vertically discretised 
    z_cav_vec = (z_wall_discr_vec(1:end-1) + z_wall_discr_vec(2:end))./2;
       
%% Building up coordinates of the stacks of cavities around the wall 
%Defining the positions of the virtual cavities along the perimeter of the
%wall
%The cavities will be located at the midpoint of each "squared section". 
% Perimeter discretization (symmetric, midpointed)
    dx = delta_xyperimeter_cavities;

% number of segments per wall (rounded so total length stays exact)
    n1 = max(1,round((2*L_x/dx)));
    n2 = max(1,round((2*L_y/dx)));
    dx1 = (2*L_x/n1);
    dx2 = (2*L_y/n2);

% --- PERIMETER NODES (exactly on the rectangle, including corners) ---
% Wall 1 (bottom): y = -L_y, x from -L_x to +L_x
    xnode_w1 = linspace(-L_x,  L_x,  n1+1);    ynode_w1 = -L_y*ones(1,n1+1);

% Wall 2 (right):  x = +L_x, y from -L_y to +L_y
    ynode_w2 = linspace(-L_y,  L_y,  n2+1);    xnode_w2 =  L_x*ones(1,n2+1);

% Wall 3 (top):    y = +L_y, x from +L_x to -L_x (reverse to keep traversal)
    xnode_w3 = linspace( L_x, -L_x, n1+1);     ynode_w3 =  L_y*ones(1,n1+1);

% Wall 4 (left):   x = -L_x, y from +L_y to -L_y
    ynode_w4 = linspace( L_y, -L_y, n2+1);     xnode_w4 = -L_x*ones(1,n2+1);

% Concatenate once per corner (avoid duplicates):
    Xnode_vec = [ xnode_w1, xnode_w2(2:end), xnode_w3(2:end), xnode_w4(2:end) ];
    Ynode_vec = [ ynode_w1, ynode_w2(2:end), ynode_w3(2:end), ynode_w4(2:end) ];

% Contiguous node index ranges per wall
    idx_node_w1 = 1 : (n1+1);
    idx_node_w2 = (n1+1) : (n1+1) + n2;
    idx_node_w3 = idx_node_w2(end) + (1:n1);
    idx_node_w4 = idx_node_w3(end) + (1:n2);

% --- CAVITY MIDPOINTS (between consecutive nodes along each wall) ---
% Wall 1 mids
    xmid_w1 = 0.5*(xnode_w1(1:end-1) + xnode_w1(2:end));
    ymid_w1 = -L_y*ones(1,n1);

% Wall 2 mids
    ymid_w2 = 0.5*(ynode_w2(1:end-1) + ynode_w2(2:end));
    xmid_w2 =  L_x*ones(1,n2);

% Wall 3 mids  (nodes reversed → mids follow traversal)
    xmid_w3 = 0.5*(xnode_w3(1:end-1) + xnode_w3(2:end));
    ymid_w3 =  L_y*ones(1,n1);

% Wall 4 mids
ymid_w4 = 0.5*(ynode_w4(1:end-1) + ynode_w4(2:end));
xmid_w4 = -L_x*ones(1,n2);

% Concatenate in the same order (w1 → w2 → w3 → w4)
Xc_vec = [xmid_w1, xmid_w2, xmid_w3, xmid_w4];
Yc_vec = [ymid_w1, ymid_w2, ymid_w3, ymid_w4];

% Contiguous cavity index ranges per wall (for your loops)
index_wall1 = 1 : n1;
index_wall2 = (n1+1) : (n1+n2);
index_wall3 = (n1+n2+1) : (2*n1+n2);
index_wall4 = (2*n1+n2+1) : (2*n1+2*n2);

% Wall lengths for normalization/labels
Lwall_1 = 2*L_x;  Lwall_3 = 2*L_x;
Lwall_2 = 2*L_y;  Lwall_4 = 2*L_y;

% --- Depth × along-wall grids ---
% Use MIDPOINT coordinates (xmid/ymid) for evaluating deflections and volumes
[Ymidwall_wall_1, Z_wall_discr_1] = meshgrid(xmid_w1, z_wall_discr_vec);
[Ymidwall_wall_2, Z_wall_discr_2] = meshgrid(ymid_w2, z_wall_discr_vec);
[Ymidwall_wall_3, Z_wall_discr_3] = meshgrid(xmid_w3, z_wall_discr_vec);
[Ymidwall_wall_4, Z_wall_discr_4] = meshgrid(ymid_w4, z_wall_discr_vec);

[Ymidwall_cav_1,  Z_cav_discr_1]  = meshgrid(xmid_w1, z_cav_vec);
[Ymidwall_cav_2,  Z_cav_discr_2]  = meshgrid(ymid_w2, z_cav_vec);
[Ymidwall_cav_3,  Z_cav_discr_3]  = meshgrid(xmid_w3, z_cav_vec);
[Ymidwall_cav_4,  Z_cav_discr_4]  = meshgrid(ymid_w4, z_cav_vec);

% Counts (driven by what we actually built)
num_cav_in_stack = numel(z_cav_vec);        % vertical bins
num_stacks       = numel(Xc_vec);           % along-wall bins (midpoints)

%% Plot results Mesh and Discretization Location
%Plots help identify the Wall Locations
if make_plots_checks == true

    figure(99)
    subplot(1,2,1); hold on; axis equal;
    plot([Xnode_vec, Xnode_vec(1)], [Ynode_vec, Ynode_vec(1)], 'b-'); % Rectangle-like shape
    plot(Xc_vec, Yc_vec, 'bo', 'LineWidth', 2, 'MarkerSize', 4,'displayname','Cavity centers'); 
    plot(Xc_vec(index_wall1), Yc_vec(index_wall1), 'cx', 'LineWidth', 2, 'MarkerSize', 8,'displayname','wall 1'); % Rectangle-like shape
    plot(Xc_vec(index_wall2), Yc_vec(index_wall2), 'y^', 'LineWidth', 2, 'MarkerSize', 8,'displayname','wall 2'); % Rectangle-like shape
    plot(Xc_vec(index_wall3), Yc_vec(index_wall3), 'mv', 'LineWidth', 2, 'MarkerSize', 8,'displayname','wall 3'); % Rectangle-like shape
    plot(Xc_vec(index_wall4), Yc_vec(index_wall4), 'k*', 'LineWidth', 2, 'MarkerSize', 8,'displayname','wall 4'); % Rectangle-like shape
    plot([-L_x L_x L_x -L_x -L_x], [-L_y -L_y L_y L_y -L_y], 'r--', 'LineWidth', 1.5, 'MarkerSize', 8,'displayname','Mesh'); % Rectangle boundary
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k','displayname','origin'); % Origin
    xlabel('X-axis'); ylabel('Y-axis');
    title('Meshing and labelling of walls');
    legend;
    grid on;
        
    subplot(1,2,2); hold on; axis equal;
    plot([Xnode_vec, Xnode_vec(1)], [Ynode_vec, Ynode_vec(1)], 'b-','displayname','Wall'); % Rectangle-like shape
    plot(Xc_vec, Yc_vec, 'bo', 'LineWidth', 2, 'MarkerSize', 4,'displayname','Cavity Centers'); % CAvity centers
    plot([-L_x L_x L_x -L_x -L_x], [-L_y -L_y L_y L_y -L_y], 'r--', 'LineWidth', 1.5, 'MarkerSize', 8,'displayname','Mesh'); % Rectangle boundary
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k','displayname','origin'); % Origin
    xlabel('X-axis'); ylabel('Y-axis');
    title('Meshing and labelling of walls');
    legend;
    grid on;

end

%% ===== Shape selection (readable handles, no giant switch) =====
% Along-wall coordinate for each wall (centered around 0 by construction)
s_w1 = xmid_w1;           % walls 1 & 3 vary in x
s_w2 = ymid_w2;           % walls 2 & 4 vary in y
s_w3 = xmid_w3;
s_w4 = ymid_w4;


if ismember(switch_shape, [3, 30, 31])
    delta_function = @(z) depth_parabolic(z, Hw);                 % −6/H*z^2 + 6z
elseif ismember(switch_shape, [5, 50, 51])
    delta_function = @(z) depth_dmm(z, Hw, C1, C2, C3);         % C1/C2/C3 combo
else
    error('switch_shape=%d not recognized.', switch_shape);
end

% --- pick the longitudinal reduction model ---
if ismember(switch_shape, [3, 5])
    long_function = @(s,Lw,Hw,HeHr) long_none(s);                 % no degradation
elseif ismember(switch_shape, [30, 50])
    long_function = @(s,Lw,Hw,HeHr) long_mu(s, Lw, Hw, HeHr);     % Mu & Huang
elseif ismember(switch_shape, [31, 51])
    long_function = @(s,Lw,Hw,HeHr) long_roboski(s, Lw, Hw, HeHr);% Roboski & Finno
else
    error('switch_shape=%d not mapped to a longitudinal model.', switch_shape);
end

% ----- Precompute 1D shapes (once) -----
% Delta shape vectors (nodes along z)
delta_wall = delta_function(z_wall_discr_vec(:));     % Nd_w × 1
delta_cav = delta_function(z_cav_vec(:));            % Nd_c × 1

% Evaluated along-wall reductions (once per wall)
R1 = long_function(s_w1, Lwall_1, Hw, He_Hwratio);   % 1 × Ns1
R2 = long_function(s_w2, Lwall_2, Hw, He_Hwratio);   % 1 × Ns2
R3 = long_function(s_w3, Lwall_3, Hw, He_Hwratio);
R4 = long_function(s_w4, Lwall_4, Hw, He_Hwratio);

% ----- Build 2D fields via outer products -----
% (depth × along-wall) with each wall’s local beta_CCS_wall_i
DELTA_wall_1 = beta_CCS_wall_1 .* (delta_wall * R1);
DELTA_wall_2 = beta_CCS_wall_2 .* (delta_wall * R2);
DELTA_wall_3 = beta_CCS_wall_3 .* (delta_wall * R3);
DELTA_wall_4 = beta_CCS_wall_4 .* (delta_wall * R4);

DELTA_cav_1  = beta_CCS_wall_1 .* (delta_cav * R1);
DELTA_cav_2  = beta_CCS_wall_2 .* (delta_cav * R2);
DELTA_cav_3  = beta_CCS_wall_3 .* (delta_cav * R3);
DELTA_cav_4  = beta_CCS_wall_4 .* (delta_cav * R4);

%%===== Integrate and compute beta/βCSS =====
VLW_runmeter_wall_1 = trapz(z_wall_discr_vec, DELTA_wall_1, 1);
VLW_runmeter_wall_2 = trapz(z_wall_discr_vec, DELTA_wall_2, 1);
VLW_runmeter_wall_3 = trapz(z_wall_discr_vec, DELTA_wall_3, 1);
VLW_runmeter_wall_4 = trapz(z_wall_discr_vec, DELTA_wall_4, 1);

beta_runmeter_wall_1 = VLW_runmeter_wall_1 / Hw^2;
beta_runmeter_wall_2 = VLW_runmeter_wall_2 / Hw^2;
beta_runmeter_wall_3 = VLW_runmeter_wall_3 / Hw^2;
beta_runmeter_wall_4 = VLW_runmeter_wall_4 / Hw^2;

beta_betaCSS_wall_1 = VLW_runmeter_wall_1 / Hw^2 / beta_CCS_wall_1;
beta_betaCSS_wall_2 = VLW_runmeter_wall_2 / Hw^2 / beta_CCS_wall_2;
beta_betaCSS_wall_3 = VLW_runmeter_wall_3 / Hw^2 / beta_CCS_wall_3;
beta_betaCSS_wall_4 = VLW_runmeter_wall_4 / Hw^2 / beta_CCS_wall_4;

%%
if make_plots_checks == true
    figure(98) % Plot the distribution of the cavities of each wall. 
    subplot(2,2,1)
    plot3(Ymidwall_wall_1/Hw,Z_wall_discr_1/Hw,DELTA_wall_1,'bo');
    xlabel('y_{wall}/Hw (-)') ; ylabel('z/H_w (-)'); zlabel('\delta_{wall} (-)') ;title('Wall 1')
    subplot(2,2,2)
    plot3(Ymidwall_wall_2/Hw,Z_wall_discr_2/Hw,DELTA_wall_2,'bo');
    xlabel('y_{wall}/Hw (-)') ; ylabel('z/H_w (-)'); zlabel('\delta_{wall} (-)') ;title('Wall 2')
    subplot(2,2,3)
    plot3(Ymidwall_wall_3/Hw,Z_wall_discr_3/Hw,DELTA_wall_3,'bo');
    xlabel('y_{wall}/Hw (-)') ; ylabel('z/H_w (-)'); zlabel('\delta_{wall} (-)') ;title('Wall 3')
    subplot(2,2,4)
    plot3(Ymidwall_wall_4/Hw,Z_wall_discr_4/Hw,DELTA_wall_4,'bo');
    xlabel('y_{wall}/Hw (-)') ; ylabel('z/H_w (-)'); zlabel('\delta_{wall} (-)') ;title('Wall 4')
    
    figure(97) % Plot the distribution of any of beta along the transverse location. Beta / Beta_css = 1 if no change from the central location. 
    subplot(2,2,1)
    plot(Ymidwall_wall_1(1,:)/Hw,beta_betaCSS_wall_1,'b-');
    xlabel('y_{wall}/Hw (-)') ; ylabel('\beta/\beta_{CSS} (-)') ;title('Wall 1')
    subplot(2,2,2)
    plot(Ymidwall_wall_2(1,:)/Hw,beta_betaCSS_wall_2,'b-');
    xlabel('y_{wall}/Hw (-)') ; ylabel('\beta/\beta_{CSS} (-)') ;title('Wall 2')
    subplot(2,2,3)
    plot(Ymidwall_wall_3(1,:)/Hw,beta_betaCSS_wall_3,'b-');
    xlabel('y_{wall}/Hw (-)') ; ylabel('\beta/\beta_{CSS} (-)') ;title('Wall 3')
    subplot(2,2,4)
    plot(Ymidwall_wall_4(1,:)/Hw,beta_betaCSS_wall_4,'b-');
    xlabel('y_{wall}/Hw (-)') ; ylabel('\beta/\beta_{CSS} (-)') ;title('Wall 4')
end
%% Locations and Depth and Deflections 
% at the wall deflection discretivaton points

Ymidwall_wall_all=[Ymidwall_wall_1,Ymidwall_wall_2,Ymidwall_wall_3,Ymidwall_wall_4]; %Lateral
Z_wall_discr_all=[Z_wall_discr_1,Z_wall_discr_2,Z_wall_discr_3,Z_wall_discr_4];     %Depths
Delta_wall_all=[DELTA_wall_1,DELTA_wall_2,DELTA_wall_3,DELTA_wall_4]; %Wall deflection magnitude

% At the cavity positions 
Ymidwall_cav_all=[Ymidwall_cav_1,Ymidwall_cav_2,Ymidwall_cav_3,Ymidwall_cav_4];
Z_cav_discr_all=[Z_cav_discr_1,Z_cav_discr_2,Z_cav_discr_3,Z_cav_discr_4];     
Delta_cav_all=[DELTA_cav_1,DELTA_cav_2,DELTA_cav_3,DELTA_cav_4];
      

%% ================= Displacement field via equivalent spherical cavities =================
% Initialisation displacement vectors 
Sx_SPn=zeros(size(Xsection)); %no symmetry, displacements from all cavity sources directly
Sy_SPn=zeros(size(Xsection));
Sz_SPn=zeros(size(Xsection));

Sx_SPy=zeros(size(Xsection)); %With symmetry plane, displacement modified to model symmetrical soil behavioy ron either side of the wall
Sy_SPy=zeros(size(Xsection));
Sz_SPy=zeros(size(Xsection));

% === Symmetry factor (scalar) ===
switch switch_solution_type
    case 1, m_symmetry = 2;   % single wall mirrored
    case 2, m_symmetry = 1;   % 4 walls analytical
    case 3, m_symmetry = 2;   % 4 walls semi-analytical (taper applied later)
    otherwise, error('Unknown switch_solution_type=%d', switch_solution_type);
end

% === Segment length per stack Δs ===
% dx1, dx2, index_wall1..4, num_cav_in_stac
seg_len = [ dx1*ones(1, numel(index_wall1)), ...
            dx2*ones(1, numel(index_wall2)), ...
            dx1*ones(1, numel(index_wall3)), ...
            dx2*ones(1, numel(index_wall4)) ];

% -- precompute equivalent cavity radii a_ij for all depths/stacks --
% ΔV_ij = Δ_cav_ij * Δz * Δs_j;  a_ij = [(3/(4π)) * m_symmetry * ΔV_ij]^(1/3)
index_wall_all  = {index_wall1, index_wall2, index_wall3, index_wall4};
DeltaS          = repmat(seg_len, num_cav_in_stack, 1);                     % [Nd_c × Ns]
vol_elem        = m_symmetry * Delta_cav_all .* (delta_z_cavities .* DeltaS);
vol_elem        = max(vol_elem, 0);                                         % clamp tiny negatives
a_cavity_matrix = ((3/(4*pi)) * vol_elem) .^ (1/3);                         % [Nd_c × Ns]



% -- quick volume consistency check --
check = [index_wall1(ceil(end/2)), index_wall2(ceil(end/2)), index_wall3(ceil(end/2)), index_wall4(ceil(end/2))];
for k = 1:numel(check)
    j = check(k);
    V_from_a     = sum((4*pi/3) * a_cavity_matrix(:,j).^3) / m_symmetry;
    V_from_delta = sum(Delta_cav_all(:,j)) * delta_z_cavities * seg_len(j);
    fprintf('[Stack %d] Vol(a)=%.3e  Vol(delta)=%.3e  rel.err=%.2e\n', j, V_from_a, V_from_delta, (V_from_a-V_from_delta)/max(1e-12,V_from_delta));
end

% -- precompute symmetry taper masks once per wall (only used if solution_type==3) --
F1 = ones(size(Xsection)); F2 = F1; F3 = F1; F4 = F1;
if switch_solution_type == 3
    % W1 (bottom, Y=-Ly): 1 at Y=-Ly → 0 at Y=+Ly
    mask = (Ysection >= -L_y);
    F1(mask) = max(0, 0.5 - 0.5*(Ysection(mask)/L_y));

    % W2 (right, X=+Lx):  1 at X=+Lx → 0 at X=-Lx
    mask = (Xsection <= +L_x);
    F2(mask) = max(0, 0.5 + 0.5*(Xsection(mask)/L_x));

    % W3 (top, Y=+Ly):    1 at Y=+Ly → 0 at Y=-Ly
    mask = (Ysection <= +L_y);
    F3(mask) = max(0, 0.5 + 0.5*(Ysection(mask)/L_y));

    % W4 (left, X=-Lx):   1 at X=-Lx → 0 at X=+Lx
    mask = (Xsection >= -L_x);
    F4(mask) = max(0, 0.5 - 0.5*(Xsection(mask)/L_x));
end
Fwall = {F1, F2, F3, F4};

% ================= accumulate displacement contributions =================
for w = 1:4
    idx_stacks = index_wall_all{w};   % stacks belonging to wall w
    Fw         = Fwall{w};            % taper mask for wall w (all ones if not used)

    for jj = idx_stacks
        % cavity stack center in plan
        deltaXl = Xc_vec(jj);
        deltaYl = Yc_vec(jj);

        % grid translated to this stack's local coordinates
        Xloc = Xsection - deltaXl;
        Yloc = Ysection - deltaYl;

        % sum depth contributions for this stack
        for ii = 1:num_cav_in_stack
            a_cavity = a_cavity_matrix(ii, jj);
            if a_cavity <= 0, continue; end

            z_cav_ii = Z_cav_discr_all(ii, jj);

            % incremental displacement from this equivalent cavity
            [ux, uy, uz] = Eq_shaft_3d_AF_v2(Xloc, Yloc, Zsection, z_cav_ii, nu, a_cavity);

            % accumulate
            Sx_SPn = Sx_SPn + ux;     Sy_SPn = Sy_SPn + uy;     Sz_SPn = Sz_SPn + uz;        % no taper
            Sx_SPy = Sx_SPy + Fw.*ux; Sy_SPy = Sy_SPy + Fw.*uy; Sz_SPy = Sz_SPy + Fw.*uz;    % with taper
        end
    end
end

% keep unmasked copies (if needed later)
Sx_SPy_full = Sx_SPy;  
Sy_SPy_full = Sy_SPy; 
Sz_SPy_full = Sz_SPy;

% zero inside the rectangle (only keep greenfield outside the box)
index_inside_rect  = (Xsection >= -L_x) & (Xsection <= L_x) & (Ysection >= -L_y) & (Ysection <= L_y);
index_outside_rect = ~index_inside_rect;

Sx_SPn = Sx_SPn .* index_outside_rect;   
Sy_SPn = Sy_SPn .* index_outside_rect;   
Sz_SPn = Sz_SPn .* index_outside_rect;

Sx_SPy = Sx_SPy .* index_outside_rect;   
Sy_SPy = Sy_SPy .* index_outside_rect;   
Sz_SPy = Sz_SPy .* index_outside_rect;


%% Extract Displacements ux,uy,uz at desired nodes or location
if switch_outputlocation == 3
    [~, iy0] = min(abs(y_vect - 0));   % row index of y≈0
    
    % === 3) Limit x to [L_x, 5*Hw] and grab that slice ===
    x_mask = (x_vect >= L_x) & (x_vect <= 5*Hw);
    x_line = x_vect(x_mask);
    
    % no-taper field along y=0
    ux_spn = Sx_SPn(iy0, x_mask).';    % column vectors for convenience
    uy_spn = Sy_SPn(iy0, x_mask).';
    uz_spn = Sz_SPn(iy0, x_mask).';
    
    % with-taper field along y=0
    ux_spy = Sx_SPy(iy0, x_mask).';
    uy_spy = Sy_SPy(iy0, x_mask).';
    uz_spy = Sz_SPy(iy0, x_mask).';
end

%% ---------------- Pack outputs (minimal) ----------------
out.Hw = Hw;
out.nu = nu;
out.L_x = L_x;
out.L_y = L_y;
out.switch_shape = switch_shape;

% grids
out.X = Xsection; 
out.Y = Ysection;
out.Z = Zsection;

% full fields (no taper and with taper)
out.Ux = Sx_SPn;
out.Uy = Sy_SPn;
out.Uz = Sz_SPn;
out.Ux_taper = Sx_SPy;
out.Uy_taper = Sy_SPy;
out.Uz_taper = Sz_SPy;

if switch_outputlocation == 3 
    % Store centerline at y≈0 (x from L_x to 5*Hw) results
    out.x_line = x_line;
    out.ux_line = ux_spn; out.uy_line = uy_spn; out.uz_line = uz_spn;
    if exist('ux_spy','var')
        out.ux_line_taper = ux_spy; 
        out.uy_line_taper = uy_spy; 
        out.uz_line_taper = uz_spy; 
    end
end 

out.meta = struct('created', datestr(now), 'mode', '3D');


if make_plots_checks == true
    if switch_outputlocation == 3

    figure(10); 
    subplot(3,1,1)
    plot(x_line./Hw, ux_spn*1000, 'k-', x_line./Hw, ux_spy*1000, 'r--','LineWidth',1.2)
    grid on; ylim([-30 30]); ylabel('u_x')
    legend('no taper (SPn)','with taper (SPy)','Location','best')
    title('y = 0 centerline')
    
    subplot(3,1,2)
    plot(x_line./Hw, uy_spn*1000, 'k-', x_line./Hw, uy_spy*1000, 'r--','LineWidth',1.2)
    grid on; ylim([-30 30]); ylabel('u_y')
    
    subplot(3,1,3)
    plot(x_line./Hw, uz_spn*1000, 'k-', x_line./Hw, uz_spy*1000, 'r--','LineWidth',1.2)
    grid on; ylim([-30 30]); ylabel('u_z'); xlabel('x/H_w (-)')
    end 
end

%% Define rectangle vertices for fill
rect_x = [-L_x L_x L_x -L_x];
rect_y = [-L_y -L_y L_y L_y];

if make_plots_results == true

    if switch_outputlocation==1
        figure(2)
        subplot(1,3,1)
        contourf(Xsection/Hw,Zsection/Hw,Sx_SPn/(beta_CCS_avg*Hw),[-1:0.1:1])
        pbaspect([1 1 1])
        colorbar; 
        caxis([-1 1 ]);
        set ( gca, 'ydir', 'reverse' )
        xlabel('x/Hw (-)') ; ylabel('z/Hw (-)') ;title('Ux_{SPn}/(\beta_{CCS}*Hw) (%)')
        
        subplot(1,3,2)
        contourf(Xsection/Hw,Zsection/Hw,Sz_SPn/(beta_CCS_avg*Hw),[-1:0.1:1])
        pbaspect([1 1 1])
        colorbar; 
        caxis([-1 1 ]);
        set ( gca, 'ydir', 'reverse' )
        xlabel('x/Hw (-)') ; ylabel('z/Hw (-)') ;title('Uz_{SPn}/(\beta_{CCS}*Hw) (%)')
        
        subplot(1,3,3)
        contourf(Xsection/Hw,Zsection/Hw,Sy_SPn/(beta_CCS_avg*Hw),[-1:0.1:1])
        pbaspect([1 1 1])
        colorbar;
        caxis([-1 1 ]);
        set ( gca, 'ydir', 'reverse' )
        xlabel('x/Hw (-)') ; ylabel('z/Hw (-)') ;title('Uy_{SPn}/(\beta_{CCS}*Hw) (%)')
        
       
        figure(20)
        subplot(1,3,1)
        contourf(Xsection/Hw,Zsection/Hw,Sx_SPy/(beta_CCS_avg*Hw),[-1:0.1:1])
        pbaspect([1 1 1])
        colorbar; 
        caxis([-1 1 ]);
        set ( gca, 'ydir', 'reverse' )
        xlabel('x/Hw (-)') ; ylabel('z/Hw (-)') ;title('Ux_{SPy}/(\beta_{CCS}*Hw) (%)')
        
        subplot(1,3,2)
        contourf(Xsection/Hw,Zsection/Hw,Sz_SPy/(beta_CCS_avg*Hw),[-1:0.1:1])
        pbaspect([1 1 1])
        colorbar; 
        caxis([-1 1 ]);
        set ( gca, 'ydir', 'reverse' )
        xlabel('x/Hw (-)') ; ylabel('z/Hw (-)') ;title('Uz_{SPy}/(\beta_{CCS}*Hw) (%)')
        
        subplot(1,3,3)
        contourf(Xsection/Hw,Zsection/Hw,Sy_SPy/(beta_CCS_avg*Hw),[-1:0.1:1])
        pbaspect([1 1 1])
        colorbar;
        caxis([-1 1 ]);
        set ( gca, 'ydir', 'reverse' )
        xlabel('x/Hw (-)') ; ylabel('z/Hw (-)') ;title('Uy_{SPy}/(\beta_{CCS}*Hw) (%)')
        
    end


    if switch_outputlocation==3
        
    figure(1); clf;
    screen_size = get(0, 'Screensize'); figure_width = screen_size(3) * 0.8; figure_height = screen_size(4) * 0.8; set(gcf, 'Position', [screen_size(1), screen_size(2), figure_width, figure_height]);
    subplot(2,3,1); hold on;
    contourf(Xsection/Hw,Ysection/Hw,Sx_SPn/(beta_CCS_avg*Hw),[-1:0.1:1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % Filled rectangle
    colorbar; 
    caxis([-1 1 ]);
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Ux_{SPn}/(\beta_{CCS}*Hw) (%)')
    set(gca,'DataAspectRatio',[1 1 1])
    
    subplot(2,3,2); hold on;
    contourf(Xsection/Hw,Ysection/Hw,Sz_SPn/(beta_CCS_avg*Hw),[-1:0.1:1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % Filled rectangle
    colorbar; 
    caxis([-1 1 ]);
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Uz_{SPn}/(\beta_{CCS}*Hw) (%)')
    set(gca,'DataAspectRatio',[1 1 1])
    
    subplot(2,2,3); hold on;
    contourf(Xsection/Hw,Ysection/Hw,Sy_SPn/(beta_CCS_avg*Hw),[-1:0.1:1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % Filled rectangle
    colorbar; 
    caxis([-1 1 ]);
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Uy_{SPn}/(\beta_{CCS}*Hw) (%)')
    set(gca,'DataAspectRatio',[1 1 1])
    
    
    subplot(2,3,3); hold on;
    contourf(Xsection/Hw,Ysection/Hw,Sy_SPn/(beta_CCS_avg*Hw),[-1:0.1:1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % Filled rectangle
    colorbar; 
    caxis([-1 1 ]);
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Uy_{SPn}/(\beta_{CCS}*Hw) (%)')
    set(gca,'DataAspectRatio',[1 1 1])
    
    subplot(2,3,4);hold on;
    contourf(Xsection/Hw,Ysection/Hw,(Sy_SPn.^2+Sx_SPn.^2).^0.5./((beta_CCS_avg*Hw)),[-1:0.1:1])
    colorbar; 
    caxis([-1 1]);
    % quiver(Xsection(1:10:end,1:10:end)/Hw,Ysection(1:10:end,1:10:end)/Hw,Sx_SPn(1:10:end,1:10:end),Sy_SPn(1:10:end,1:10:end), 'k', 'LineWidth', 1.5)
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('u_h_{SPn}/(\beta_{CCS} Hw)=sqrt(ux_{SPn}^2+uy_{SPn}^2)/(\beta_{CCS} Hw)')
    set(gca,'DataAspectRatio',[1 1 1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    
    subplot(2,3,5);hold on;
    contourf(Xsection/Hw,Ysection/Hw,(Sy_SPn.^2+Sx_SPn.^2).^0.5./Sz_SPn,[0:0.5:10])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    colorbar; 
    caxis([0 5 ]);
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('u_h_{SPn}/u_z_{SPn}')
    set(gca,'DataAspectRatio',[1 1 1])
    
    
    figure(2);clf; 
    screen_size = get(0, 'Screensize'); figure_width = screen_size(3) * 0.8; figure_height = screen_size(4) * 0.8; set(gcf, 'Position', [screen_size(1), screen_size(2), figure_width, figure_height]);
    subplot(2,3,1);hold on;
    contourf(Xsection/Hw,Ysection/Hw,atand(Sy_SPn./Sx_SPn),[-90:10:90])
    colorbar; 
    caxis([-90 90 ]);
    quiver(Xsection(1:10:end,1:10:end)/Hw,Ysection(1:10:end,1:10:end)/Hw,Sx_SPn(1:10:end,1:10:end),Sy_SPn(1:10:end,1:10:end), 'k', 'LineWidth', 1.5)
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Omega_{SPn} (deg) from analytical model')
    set(gca,'DataAspectRatio',[1 1 1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    
    
    subplot(2,3,2);hold on;
    % % contourf(Xsection/Hw,Ysection/Hw,atand(Ysection./Xsection),[-90:10:90])
    % % colorbar; 
    % % caxis([-90 90 ]);
    % % quiver(Xsection(1:10:end,1:10:end)/Hw,Ysection(1:10:end,1:10:end)/Hw,Sx_SPn(1:10:end,1:10:end),Sy_SPn(1:10:end,1:10:end), 'k', 'LineWidth', 1.5)
    % % set ( gca, 'ydir', 'reverse' )
    % % xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Omega_{SPn} (deg) from empirical assumption')
    % % set(gca,'DataAspectRatio',[1 1 1])
    % % fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    contour(Xsection/Hw,Ysection/Hw,atand(Sy_SPn./Sx_SPn),    [-90:15:90],'-k',"ShowText",true,'displayname','\omega')
    [Contour1, h1] =contour(Xsection/Hw,Ysection/Hw,atand(Ysection./Xsection),[-90:15:90],'--r',"ShowText",true,'displayname','u_h towards centre');
    clabel(Contour1, h1, 'Color', 'r'); % Set label color to black
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Analysis \omega SPn')
    set(gca,'DataAspectRatio',[1 1 1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5,'displayname','station'); % ,'facecolor','white'
    legend('location','best')
    
    
    
    subplot(2,3,3);hold on;
    contour(Xsection/Hw,Ysection/Hw,atand(Sy_SPn./Sx_SPn),[-90:15:90],'-k',"ShowText",true,'displayname','\omega')
    contour(Xsection/Hw,Ysection/Hw,Sz_SPn/(beta_CCS_avg*Hw),[1, 2, 5, 10, 50, 100]/100,'--r',"ShowText",true,'displayname','Uz}/\beta/H_w')
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Analysis \omega SPn')
    set(gca,'DataAspectRatio',[1 1 1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5,'displayname','station'); % ,'facecolor','white'
    legend('location','best')
    
    
    
    subplot(2,3,4);hold on;
    contourf(Xsection/Hw,Ysection/Hw,(Sy_SPn.^2+Sx_SPn.^2).^0.5./((beta_CCS_avg*Hw)),[0:0.1:1])
    colorbar; 
    caxis([0 1]);
    quiver(Xsection(1:10:end,1:10:end)/Hw,Ysection(1:10:end,1:10:end)/Hw,Sx_SPn(1:10:end,1:10:end),Sy_SPn(1:10:end,1:10:end), 'k', 'LineWidth', 1.5)
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('u_h_{SPn}/(\beta_{CCS} Hw)=sqrt(ux_{SPn}^2+uy_{SPn}^2)/(\beta_{CCS} Hw)')
    set(gca,'DataAspectRatio',[1 1 1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    
    subplot(2,3,5);hold on;
    colorbar; 
    caxis([0 1]);
    quiver(Xsection(1:10:end,1:10:end)/Hw,Ysection(1:10:end,1:10:end)/Hw,Sx_SPn(1:10:end,1:10:end),Sy_SPn(1:10:end,1:10:end), 'k', 'LineWidth', 1.5)
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('u_{tot}_{SPn}/(\beta_{CCS} Hw)=sqrt(ux_{SPn}^2+uy_{SPn}^2+u_z_{SPn}^2)/(\beta_{CCS} Hw)')
    set(gca,'DataAspectRatio',[1 1 1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    contourf(Xsection/Hw,Ysection/Hw,(Sy_SPn.^2+Sx_SPn.^2+Sz_SPn.^2).^0.5./((beta_CCS_avg*Hw)),[-1:0.1:1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    figure(10); clf;
    screen_size = get(0, 'Screensize'); figure_width = screen_size(3) * 0.8; figure_height = screen_size(4) * 0.8; set(gcf, 'Position', [screen_size(1), screen_size(2), figure_width, figure_height]);
    subplot(2,3,1); hold on;
    contourf(Xsection/Hw,Ysection/Hw,Sx_SPy/(beta_CCS_avg*Hw),[-1:0.1:1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % Filled rectangle
    colorbar; 
    caxis([-1 1 ]);
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Ux_{SPy}/(\beta_{CCS}*Hw) (%)')
    set(gca,'DataAspectRatio',[1 1 1])
    
    subplot(2,3,2); hold on;
    contourf(Xsection/Hw,Ysection/Hw,Sz_SPy/(beta_CCS_avg*Hw),[-1:0.1:1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % Filled rectangle
    colorbar; 
    caxis([-1 1 ]);
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Uz_{SPy}/(\beta_{CCS}*Hw) (%)')
    set(gca,'DataAspectRatio',[1 1 1])
    
    subplot(2,3,3); hold on;
    contourf(Xsection/Hw,Ysection/Hw,Sy_SPy/(beta_CCS_avg*Hw),[-1:0.1:1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % Filled rectangle
    colorbar; 
    caxis([-1 1 ]);
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Uy_{SPy}/(\beta_{CCS}*Hw) (%)')
    set(gca,'DataAspectRatio',[1 1 1])
    
    
    subplot(2,3,4);hold on;
    contourf(Xsection/Hw,Ysection/Hw,(Sy_SPy.^2+Sx_SPy.^2).^0.5./((beta_CCS_avg*Hw)),[-1:0.1:1])
    colorbar; 
    caxis([-1 1]);
    % quiver(Xsection(1:10:end,1:10:end)/Hw,Ysection(1:10:end,1:10:end)/Hw,Sx_SPn(1:10:end,1:10:end),Sy_SPn(1:10:end,1:10:end), 'k', 'LineWidth', 1.5)
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('u_h_{SPy}/(\beta_{CCS} Hw)=sqrt(ux_{SPy}^2+uy_{SPn}^2)/(\beta_{CCS} Hw)')
    set(gca,'DataAspectRatio',[1 1 1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    
    subplot(2,3,5);hold on;
    contourf(Xsection/Hw,Ysection/Hw,(Sy_SPy.^2+Sx_SPy.^2).^0.5./Sz_SPy,[0:0.5:10])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    colorbar; 
    caxis([0 5 ]);
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('u_h_{SPy}/u_z_{SPy}')
    set(gca,'DataAspectRatio',[1 1 1])
    
    
    figure(20);clf; 
    screen_size = get(0, 'Screensize'); figure_width = screen_size(3) * 0.8; figure_height = screen_size(4) * 0.8; set(gcf, 'Position', [screen_size(1), screen_size(2), figure_width, figure_height]);
    subplot(2,3,1);hold on;
    contourf(Xsection/Hw,Ysection/Hw,atand(Sy_SPy./Sx_SPy),[-90:10:90])
    colorbar; 
    caxis([-90 90 ]);
    quiver(Xsection(1:10:end,1:10:end)/Hw,Ysection(1:10:end,1:10:end)/Hw,Sx_SPy(1:10:end,1:10:end),Sy_SPy(1:10:end,1:10:end), 'k', 'LineWidth', 1.5)
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Omega_{SPy} (deg) from analytical model')
    set(gca,'DataAspectRatio',[1 1 1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    
    subplot(2,3,2);hold on;
    % contourf(Xsection/Hw,Ysection/Hw,atand(Ysection./Xsection),[-90:10:90])
    % fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    % colorbar; 
    % caxis([-90 90 ]);
    % quiver(Xsection(1:10:end,1:10:end)/Hw,Ysection(1:10:end,1:10:end)/Hw,Sx_SPy(1:10:end,1:10:end),Sy_SPy(1:10:end,1:10:end), 'k', 'LineWidth', 1.5)
    % set ( gca, 'ydir', 'reverse' )
    % xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Omega_{SPy} (deg) from empirical assumption')
    % set(gca,'DataAspectRatio',[1 1 1])
    contour(Xsection/Hw,Ysection/Hw,atand(Sy_SPy./Sx_SPy),    [-90:15:90],'-k',"ShowText",true,'displayname','\omega')
    [Contour1, h1] =contour(Xsection/Hw,Ysection/Hw,atand(Ysection./Xsection),[-90:15:90],'--r',"ShowText",true,'displayname','u_h towards centre');
    clabel(Contour1, h1, 'Color', 'r'); % Set label color to black
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Analysis \omega SPy')
    set(gca,'DataAspectRatio',[1 1 1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5,'displayname','station'); % ,'facecolor','white'
    legend('location','best')
    
    
    
    subplot(2,3,3);hold on;
    contour(Xsection/Hw,Ysection/Hw,atand(Sy_SPy./Sx_SPy),[-90:15:90],'-k',"ShowText",true,'displayname','\omega')
    contour(Xsection/Hw,Ysection/Hw,Sz_SPy/(beta_CCS_avg*Hw),[1, 2, 5, 10, 50, 100]/100,'--r',"ShowText",true,'displayname','Uz}/\beta/H_w')
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('Analysis \omega SPy')
    set(gca,'DataAspectRatio',[1 1 1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5,'displayname','station'); % ,'facecolor','white'
    legend('location','best')
    
    
    
    subplot(2,3,4);hold on;
    contourf(Xsection/Hw,Ysection/Hw,(Sy_SPy.^2+Sx_SPy.^2).^0.5./((beta_CCS_avg*Hw)),[0:0.1:1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    colorbar; 
    caxis([0 1]);
    quiver(Xsection(1:10:end,1:10:end)/Hw,Ysection(1:10:end,1:10:end)/Hw,Sx_SPy(1:10:end,1:10:end),Sy_SPy(1:10:end,1:10:end), 'k', 'LineWidth', 1.5)
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('u_h_{SPy}/(\beta_{CCS} Hw)=sqrt(ux^2+uy^2)/(\beta_{CCS} Hw)')
    set(gca,'DataAspectRatio',[1 1 1])
    
    subplot(2,3,5);hold on;
    contourf(Xsection/Hw,Ysection/Hw,(Sy_SPy.^2+Sx_SPy.^2+Sz_SPy.^2).^0.5./((beta_CCS_avg*Hw)),[0:0.1:1])
    fill(rect_x/Hw, rect_y/Hw, [0.8 0.8 0.8],'facecolor','white', 'EdgeColor', 'r', 'LineWidth', 1.5); % ,'facecolor','white'
    colorbar; 
    caxis([0 1]);
    quiver(Xsection(1:10:end,1:10:end)/Hw,Ysection(1:10:end,1:10:end)/Hw,Sx_SPy(1:10:end,1:10:end),Sy_SPy(1:10:end,1:10:end), 'k', 'LineWidth', 1.5)
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('y/Hw (-)') ;title('u_{tot}_{SPy}/(\beta_{CCS} Hw)=sqrt(ux_{SPy}^2+uy_{SPy}^2+u_z_{SPy}^2)/(\beta_{CCS} Hw)')
    set(gca,'DataAspectRatio',[1 1 1])
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%
    
    
    figure(100); clf;
    
    [a,index_centralsection_AAalongx]=min(abs(Ysection(:,1)));
    [a,index_centralsection_BBalongy]=min(abs(Xsection(1,:)));
    
    
    subplot(2,3,1); hold on;
    plot(Xsection(index_centralsection_BBalongy,:)/Hw,Sx_SPn(index_centralsection_BBalongy,:)/(beta_CCS_avg*Hw),'k-','displayname','SPn')
    plot(Xsection(index_centralsection_BBalongy,:)/Hw,Sx_SPy(index_centralsection_BBalongy,:)/(beta_CCS_avg*Hw),'r-','displayname','SPy')
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('Ux_{SPn}/(\beta_{CCS}*Hw) (%)')
    legend
    if switch_validation_plot==1
    plot(x_H_2d+L_x/Hw,Sx0_Hbeta,'m--','displayname','2D','linewidth',3)
    xlim([L_x/Hw,max(Xsection(index_centralsection_BBalongy,:)/Hw)])
    end
    
    subplot(2,3,2); hold on;
    plot(Xsection(index_centralsection_BBalongy,:)/Hw,Sz_SPn(index_centralsection_BBalongy,:)/(beta_CCS_avg*Hw),'k-','displayname','SPn')
    plot(Xsection(index_centralsection_BBalongy,:)/Hw,Sz_SPy(index_centralsection_BBalongy,:)/(beta_CCS_avg*Hw),'r-','displayname','SPy')
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('Uz_{SPn}/(\beta_{CCS}*Hw) (%)')
    legend
    if switch_validation_plot==1
    plot(x_H_2d+L_x/Hw,Sz0_Hbeta,'m--','displayname','2D','linewidth',3)  
    xlim([L_x/Hw,max(Xsection(index_centralsection_BBalongy,:)/Hw)])
    end
    
    subplot(2,3,3); hold on;
    plot(Xsection(index_centralsection_BBalongy,:)/Hw,Sy_SPn(index_centralsection_BBalongy,:)/(beta_CCS_avg*Hw),'k-','displayname','SPn')
    plot(Xsection(index_centralsection_BBalongy,:)/Hw,Sy_SPy(index_centralsection_BBalongy,:)/(beta_CCS_avg*Hw),'r-','displayname','SPy')
    set ( gca, 'ydir', 'reverse' )
    xlabel('x/Hw (-)') ; ylabel('Uy_{SPn}/(\beta_{CCS}*Hw) (%)')
    legend
    if switch_validation_plot==1
    xlim([L_x,max(Xsection(index_centralsection_BBalongy,:)/Hw)])
    end
    
    subplot(2,3,4); hold on;
    plot(Ysection(:,index_centralsection_AAalongx)/Hw,Sx_SPn(:,index_centralsection_AAalongx)/(beta_CCS_avg*Hw),'k-','displayname','SPn')
    plot(Ysection(:,index_centralsection_AAalongx)/Hw,Sx_SPy(:,index_centralsection_AAalongx)/(beta_CCS_avg*Hw),'r-','displayname','SPy')
    set ( gca, 'ydir', 'reverse' )
    xlabel('y/Hw (-)') ; ylabel('Ux_{SPn}/(\beta_{CCS}*Hw) (%)')
    legend
    
    subplot(2,3,5); hold on;
    plot(Ysection(:,index_centralsection_AAalongx)/Hw,Sz_SPn(:,index_centralsection_AAalongx)/(beta_CCS_avg*Hw),'k-','displayname','SPn')
    plot(Ysection(:,index_centralsection_AAalongx)/Hw,Sz_SPy(:,index_centralsection_AAalongx)/(beta_CCS_avg*Hw),'r-','displayname','SPy')
    set ( gca, 'ydir', 'reverse' )
    xlabel('y/Hw (-)') ; ylabel('Uz_{SPn}/(\beta_{CCS}*Hw) (%)')
    legend
    
    subplot(2,3,6); hold on;
    plot(Ysection(:,index_centralsection_AAalongx)/Hw,Sy_SPn(:,index_centralsection_AAalongx)/(beta_CCS_avg*Hw),'k-','displayname','SPn')
    plot(Ysection(:,index_centralsection_AAalongx)/Hw,Sy_SPy(:,index_centralsection_AAalongx)/(beta_CCS_avg*Hw),'r-','displayname','SPy')
    set ( gca, 'ydir', 'reverse' )
    xlabel('y/Hw (-)') ; ylabel('Uy_{SPn}/(\beta_{CCS}*Hw) (%)')
    legend
    
   
    titles = {'W1 bottom','W2 right','W3 top','W4 left'};
    figure(10000); clf; tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for w = 1:4
        nexttile; 
        F = (switch_solution_type==3) * Fwall{w} + (switch_solution_type~=3) * ones(size(Xsection));
        contourf(Xsection/Hw, Ysection/Hw, F, 20, 'LineColor','none'); axis equal tight
        hold on; fill(rect_x/Hw, rect_y/Hw, [1 1 1], 'EdgeColor','r', 'LineWidth', 1.0);
        title(titles{w}); xlabel('x/H_w'); ylabel('y/H_w');
    end
    
    
    
    end

end

%%
if make_plots_checks == true

    figure(393);clf;
    for He_L=[13.6/20, 13.6/80]
    He_Hwratio=0.5;
    
    Hw=13.6/He_Hwratio;
    He=Hw*He_Hwratio;
    
    Lwall=He/He_L;
    y_midwall=-Lwall/2:Lwall/500:Lwall/2;
   
    red_finno=  (1-0.5.*erfc( ...
                      (2.8.*((Lwall/2-abs(y_midwall))+Lwall.*(0.015+0.035.*log(He_Hwratio.*Hw./Lwall))))       ...
         ./ (0.5.*Lwall-Lwall.*(0.15+0.035.*log(He_Hwratio.*Hw/Lwall))))  ) ; 
    
        red_hu= exp(-1*pi.*(y_midwall./ (Lwall./2.*(0.069.*log(Hw./He_Hwratio./Lwall)+1.03))).^2);     
     
    % % %         tttt=0.1;
    % % %          red_wall=  (1-0.5.*erfc( ...
    % % %                   (2.8.*((Lwall/2-abs(y_midwall))+Lwall.*(0.015+0.035.*log(tttt.*Hw./Lwall))))       ...
    % % %      ./ (0.5.*Lwall-Lwall.*(0.15+0.035.*log(tttt.*Hw/Lwall))))  ) ; 
    
    subplot(1,2,1); hold on; 
    plot(y_midwall,red_finno,'-','displayname', sprintf('Finno - He/Lwall: %.2f', He_Hwratio*Hw/Lwall))
    plot(y_midwall,red_hu   ,'--','displayname', sprintf('Mu - He/Lwall: %.2f', He_Hwratio*Hw/Lwall))
    title(sprintf('Hw: %.2f m He:%.2f m', Hw,He))
    xlabel('y (m)') ; ylabel('Finno reduction factor')
    subplot(1,2,2); hold on; 
    plot(y_midwall/Lwall,red_finno,'-','displayname', sprintf('Finno - He/Lwall: %.2f', He_Hwratio*Hw/Lwall))
    plot(y_midwall/Lwall,red_hu   ,'--','displayname', sprintf('Mu - He/Lwall: %.2f', He_Hwratio*Hw/Lwall))
    title(sprintf('Hw: %.2f m He:%.2f m', Hw,He))
    xlabel('y/Lwall (-)') ; ylabel('Finno reduction factor')
    end
    % % % plot(y_midwall/Lwall,red_wall   ,'r-','displayname', 'FInno adjusted')
    legend('location','best')
end 
%% Idea, omega may be perpendicular to isolines
if switch_outputlocation == 3
    %close all
    [M_Sz_SPn,c_Sz_SPn] =contour(Xsection/Hw,Ysection/Hw,Sz_SPn/(beta_CCS_avg*Hw),[-1:0.1:1]);
    
    % Initialize storage for X, Y coordinates
    x_isoline = [];
    y_isoline = [];
    
    % Loop through contour data
    k = 1;
    while k < size(M_Sz_SPn, 2)
        level = M_Sz_SPn(1, k);    % Contour level
        num_points = M_Sz_SPn(2, k); % Number of points in this contour
        
        if level == 10  % Extract only this level  isoline
            x_isoline = [x_isoline, NaN, M_Sz_SPn(1, k+1:k+num_points)];
            y_isoline = [y_isoline, NaN, M_Sz_SPn(2, k+1:k+num_points)];
        end
    
        % Move to the next contour segment
        k = k + num_points + 1;
    end
    
    
    % Display extracted isoline points
    fprintf('Extracted %d points for the isoline at level -50\n', length(x_isoline));
    
    if make_plots_checks == true
        % Plot results
        figure(1001);hold on;
        plot(x_isoline, y_isoline, 'r', 'LineWidth', 2); % Extracted isoline in red
        for kk=1:10:length((x_isoline))
          plot([0,x_isoline(kk)],[0, y_isoline(kk)], 'k-', 'LineWidth', 2); % Extracted isoline in red  
        end
        xlabel('x/Hw (-)'); ylabel('y/Hw (-)');
        title('Extracted Isoline');
        set(gca, 'ydir', 'reverse');
        grid on;
    
    end 
end 
% ================== Local helpers ==================
%v1 has settlements positive downwarda and horizontal movements positive
%towards the cavity
%v2 has settlements positive downwarda and horizontal movements positive
%according to axes

%Pinto Solution
function [u_x, u_y, u_z]=Eq_shaft_3d_AF_v2(x,y,z,h,nu,a)  

    R1=(x.^2+y.^2+(z-h).^2).^0.5;
    R2=(x.^2+y.^2+(z+h).^2).^0.5;
    
    % single spherity 
    fx=x.*(1./R1.^3+(3-4.*nu)./R2.^3-6.*z.*(z+h)./R2.^5);
    u_x=1/3*a^3*fx;
    u_x=u_x*-1;
    
    fy=y.*(1./R1.^3+(3-4.*nu)./R2.^3-6.*z.*(z+h)./R2.^5);
    u_y=1/3*a^3*fy;
    u_y=u_y*-1;
    
    fz=-1*((z-h)./R1.^3+    2.*z./R2.^3    -(3-4.*nu).*(z+h)./R2.^3    -6.*z.*(z+h).^2./R2.^5);
    u_z=1/3*a^3*fz;

end


function delta_w = depth_parabolic(z_vec, Hw)
% Parabolic vertical shape: D(z) = -6/Hw*z^2 + 6*z  (clamped to ≥0)
    z = z_vec(:);
    delta_w = (-6./Hw .* z.^2 + 6.*z);
    delta_w = max(delta_w , 0);
end

function delta_w = depth_dmm(z_vec, Hw, C1, C2, C3)
% C1/C2/C3 vertical combo (no longitudinal degradation built-in here):
% C1: cantilever (2*Hw - 2*z), C2: parabolic, C3: kick-in (2*z)
    z = z_vec(:);
    delta_w = C1.*(2*Hw - 2*z) + C2.* (-6./Hw .* z.^2 + 6.*z) + C3.*(2*z);
    delta_w = max(delta_w, 0);
end

function reduction = long_none(s_vec)
% No longitudinal degradation
   reduction = ones(1, numel(s_vec));
end

function reduction = long_mu(s_vec, Lwall, Hw, He_Hwratio)
% Mu & Huang (2016) Gaussian reduction
% R(s) = exp( -pi * (s/W)^2 ), with W from their expression
    s = s_vec(:).';
    W = (Lwall/2) * (0.069.*log( (Hw./He_Hwratio)./Lwall ) + 1.03);
    W = max(W, 0.05*Lwall);  % safety floor
    reduction = exp( -pi .* (s./W).^2 );
end

function reduction = long_roboski(s_vec, Lwall, Hw, He_Hwratio)
% Roboski & Finno (2006) erfc-type reduction
    s = s_vec(:).';
    a = 2.8;
    c = 0.015 + 0.035*log(He_Hwratio*Hw/Lwall);
    d = 0.15  + 0.035*log(He_Hwratio*Hw/Lwall);
    num = a .* ( (Lwall/2 - abs(s)) + Lwall.*c );
    den = max( (0.5*Lwall - Lwall.*d), 1e-9 );
    reduction = 1 - 0.5 .* erfc(num ./ den);
end


function v = getf(s, name, default)
if isfield(s, name), v = s.(name);
else
    v = default; 
end
end

end

