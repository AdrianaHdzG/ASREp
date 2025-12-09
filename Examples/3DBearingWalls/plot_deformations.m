clear
close all

% ---- Load mesh + displacements ----
load('./input/kk_0811220251203T085216.mat');   % wholeNodesXYZ
ep_u = readmatrix('./output/nodal_displacement_elastoplastic.txt');

% ---- Coordinates ----
x = wholeNodesXYZ(:,1);
y = wholeNodesXYZ(:,2);
z = wholeNodesXYZ(:,3);

% ---- Displacements (to mm; remove *1000 if already in mm) ----
ux = ep_u(1:3:end) * 1000;
uy = ep_u(2:3:end) * 1000;
uz = ep_u(3:3:end) * 1000;

% ---- Marker size (adaptive with node count) ----
n = numel(x);
msize = max(6, min(60, 30000 / max(n,1)));  % small models -> bigger markers

% ---- (Optional) symmetric color limits per component for fair comparison ----
clx = [-20 1]
cly = [-20 1]
clz = [-20 1]

% ---- Plot tiles ----
figure('Name','3D nodal displacement components','Position',[100 100 1500 520])
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

% u_x
nexttile
scatter3(x, y, z, msize, ux, 'filled');
axis equal tight; grid on; box on; view(45,30)
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('u_x (mm)', 'Interpreter','none'); colormap(turbo); cb=colorbar;
caxis(clx);

% u_y
nexttile
scatter3(x, y, z, msize, uy, 'filled');
axis equal tight; grid on; box on; view(45,30)
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('u_y (mm)', 'Interpreter','none'); colormap(turbo); colorbar;
caxis(cly);

% u_z
nexttile
scatter3(x, y, z, msize, uz, 'filled');
axis equal tight; grid on; box on; view(45,30)
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('u_z (mm)', 'Interpreter','none'); colormap(turbo); colorbar;
caxis(clz);
