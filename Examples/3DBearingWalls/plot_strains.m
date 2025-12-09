clear 
close all
% Append paths to plot functions
addpath('../../Postprocessing_in_FEM')
% Load mesh
load('./input/kk_0811220251203T085216.mat');
% Load tensile strain
el_strains = readmatrix('./output/principal_strain_elastic.txt');
ep_strains = readmatrix('./output/principal_strain_elastoplastic.txt');
% Load nodal disp
el_u = readmatrix('./output/nodal_displacement_elastic.txt');
ep_u = readmatrix('./output/nodal_displacement_elastoplastic.txt');

%------ Deep Excavation Required Input
% L_x            = 9.5;     % half-length of station in x (example; use your value)
% building_offset = 11;   % distance from the wall face to the building line
% y_build_center  = 0.0;    % change if the building is not centered in y
% 
% % Shift the *entire* building (all nodes, including interface nodes):
% x_shift = (2*L_x + building_offset);
% y_shift = y_build_center;

% wholeNodesXYZ(:,1) = wholeNodesXYZ(:,1) + x_shift;



% Plot elastic strain
figure
PlotFieldonDefoMesh(wholeNodesXYZ, wholeElem2n, 100,...
    [el_u(1:3:end),el_u(2:3:end),el_u(3:3:end)],el_strains(:,1)*1000000,...
    'Principal tensile strains($\mu\varepsilon$)')
daspect([1 1 1])
%saveas(gcf, './output/strain_plot_elastic.png')
%saveas(gcf, './output/strain_plot_elastic.fig')
% Plot elastoplastic strain
figure
PlotFieldonDefoMesh(wholeNodesXYZ, wholeElem2n, 100,...
    [ep_u(1:3:end),ep_u(2:3:end),ep_u(3:3:end)],ep_strains(:,1)*1000000,...
    'Principal tensile strains($\mu\varepsilon$)')
daspect([1 1 1])
%saveas(gcf, './output/strain_plot_elastoplastic.png')
%saveas(gcf, './output/strain_plot_elastoplastic.fig')


