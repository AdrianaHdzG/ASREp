clear
close all
load('../Examples/3DBearingWalls/input/burd_mat_20250901T162943.mat')
% If using Burd 2022 or Yiu 2019 FEM simulated greenfield displacement
% Yiu_GF = load('../Examples/3DBearingWalls/input/Yiu_GF.mat');
% Yiu_GF = Yiu_GF.Yiu_GF;

wholeNodesXYZ(:,3) = wholeNodesXYZ(:,3) - min(wholeNodesXYZ(:,3));
interNodesXYZ(:,3) = 0;


%------ Deep Excavation Required Input
L_x            = 9.5;     % half-length of station in x (example; use your value)
building_offset = 11;   % distance from the wall face to the building line
y_build_center  = 0.0;    % change if the building is not centered in y

% Shift the *entire* building (all nodes, including interface nodes):
x_shift = (2*L_x + building_offset);
y_shift = y_build_center;

wholeNodesXYZ(:,1) = wholeNodesXYZ(:,1) + x_shift;
wholeNodesXYZ(:,2) = wholeNodesXYZ(:,2) + y_shift;

interNodesXYZ(:,1) = interNodesXYZ(:,1) + x_shift;
interNodesXYZ(:,2) = interNodesXYZ(:,2) + y_shift;
%------------------------------------------------

PlotMesh(wholeNodesXYZ, wholeElem2n, 0)
PlotMesh(interNodesXYZ(:, 1:2), elem2nInter, 0); daspect([1, 1, 1])

x = wholeNodesXYZ(:,1)';
y = wholeNodesXYZ(:,2)';
z = wholeNodesXYZ(:,3)';

%%
tic
%------------------------Define tunnel and building parameters--------------------
num_incr_TUNNEL = 20;
num_incr_P = 10;
width_Kx = 0.57;%Pickhaver=0.75/Yiu=0.57
width_Ky = 0.57;
tunnel_z0 = 22;%Pickhaver=10/Yiu=23-1 to model the embedment of foundation
tunnel_d = 11;%Pickhaver=5/Yiu=11
tunnel_delta = 0.3;%S_tunnelFace/S_max
tunnel_ys = -50;%location of tunnel face 
tunnel_yf = 500;%location of tunnel start point
Vltp = 1.65; %Pickhaver=1.36925/Yiu=1.65
masonry_nu = 0.2;
masonry_E = 3000*10^6;%Pickhaver=10^10/Yiu=3*10^9
masonry_G = masonry_E/2/(1+masonry_nu);
% masonry_rho =1800*9.81+10000/0.2/3;
masonry_rho = 23.75*10^3; %burd table 1
timber_E = masonry_E;
timber_nu = 0.2;
timber_G = timber_E/2/(1+timber_nu);
timber_rho = masonry_rho;
nus = 0.2;
Gs = 6000000;
Es = Gs*2*(1+nus);
mu = 0.3;
%%
NoOfElem = size(wholeElem2n,1);
eleRho = repmat([masonry_rho],1,NoOfElem);
eleE = repmat([masonry_E],1,NoOfElem);
elev = repmat([masonry_nu],1,NoOfElem);
eleG = repmat([masonry_G],1,NoOfElem);
eleRho(timberElemIndex(1):timberElemIndex(2)) = timber_rho;
eleE(timberElemIndex(1):timberElemIndex(2)) = timber_E;
elev(timberElemIndex(1):timberElemIndex(2)) = timber_nu;
eleG(timberElemIndex(1):timberElemIndex(2)) = timber_G;
[KKfoot,Q] = computeKKfootAndBodyForce_v2(x', y', z', wholeElem2n, eleE, elev, eleRho, eleG);

%----------------Calculate and plot greenfield displacement--------------------------------------------------
newNode = unique(reshape(wholeElem2n, length(wholeElem2n)*8,1));
newNodeDOF = newNode*3-3+1;
newNodeDOF = [newNodeDOF; newNode*3-3+2];
newNodeDOF = [newNodeDOF; newNode*3-3+3];
newNodeDOF = sort(newNodeDOF);
openKKfoot = KKfoot(newNodeDOF, newNodeDOF);
openQ = Q(newNodeDOF);



% [Uffx,Uffy,Uffz] = u_3D_comas_v2(Vltp,...
%         width_Kx, width_Ky, x(newNode), y(newNode), z(newNode), tunnel_delta, tunnel_z0,...
%         tunnel_ys, tunnel_yf, tunnel_d, find(z(newNode)==min(z(newNode))));

paramsGF = struct( ...
    'nu',0.499,'Hw',19,'L_x',9.5,'L_y',32, ...
    'switch_shape',5,'C1',0,'C2',1,'C3',0, ...
    'beta_CCS_wall_1',7.5e-4,'beta_CCS_wall_2',7.5e-4, ...
    'beta_CCS_wall_3',7.5e-4,'beta_CCS_wall_4',7.5e-4, ...
    'delta_z_cavities',19/19,'delta_xyperimeter_cavities',2.5, ...
    'switch_outputlocation',100, ...                % <— new mode
    'x_nodes', x(newNode)', 'y_nodes', y(newNode)', 'z_nodes', z(newNode)' );

outGF = run_greenfield_3D_DMM(paramsGF);           % run at node cloud
Uffx  = outGF.Ux(:);                               % pick no-taper components
Uffy  = outGF.Uy(:);
Uffz  = -1*outGF.Uz(:);

figure(1234); 
plot3(x(newNode), y(newNode), Uffz, '.'); grid on
xlabel x; ylabel y; zlabel 'Uffz (m)'; title('GF at building nodes (DMM)')


ox = x(newNode);
oy = y(newNode);
oz = z(newNode);

% Use Burd/Yiu GF
%Uffz(oz==min(z)) = (interp1(Yiu_GF(:,1), Yiu_GF(:,2)/(-1000), (ox(oz==min(z))).'))';
%Uffx(oz==min(z)) = (ox(oz==min(z))).'./(tunnel_z0-(oz(oz==min(z))).').*Uffz(oz==min(z)); 
%Uffy = zeros(size(Uffz));

figure; 
plot3(ox(oz==min(z)), oy(oz==min(z)), Uffz(oz==min(z)),'x')
xlabel('x'); ylabel('y'); zlabel('z (m)')
title('UFFZ')
figure;
plot3(ox(oz==min(z)), oy(oz==min(z)), Uffx(oz==min(z),1),'x')
xlabel('x'); ylabel('y'); zlabel('z (m)')
title('UFFx')
figure;
plot3(ox(oz==min(z)), oy(oz==min(z)), Uffy(oz==min(z),1),'x')
xlabel('x'); ylabel('y'); zlabel('z (m)')
title('UFFy')


bottom_nodes = find(oz==min(oz));
groundNode = bottom_nodes;
groundNodeDOF1 = bottom_nodes*3-3+1;
groundNodeDOF2 = bottom_nodes*3-3+2;
groundNodeDOF3 = bottom_nodes*3-3+3;
groundNodeDOF = [groundNodeDOF1, groundNodeDOF2, groundNodeDOF3];
groundNodeDOF = sort(groundNodeDOF);
bottom_elems = [];
for i = 1:length(wholeElem2n)
    if z(wholeElem2n(i,1)) == min(oz)
        bottom_elems = [bottom_elems;i];
    end
end
PlotMesh(wholeNodesXYZ(:,1:2), wholeElem2n(bottom_elems, 1:4), 0);
daspect([1 1 1])
%%
%----------------- trival test 1--------------------------------------
% topMove = 0.25;
% topNodes = find(oz == max(oz));
% topNodesDOF1 = topNodes*3-3+1;
% topNodesDOF2 = topNodes*3-3+2;
% topNodesDOF3 = topNodes*3-3+3;
% topNodesDOF = [topNodesDOF1, topNodesDOF2, topNodesDOF3];
% topNodesDOF = sort(topNodesDOF);
% free_nodes = find(oz~=0&oz~=max(oz));
% free_nodesDOF1 = free_nodes*3-3+1;
% free_nodesDOF2 = free_nodes*3-3+2;
% free_nodesDOF3 = free_nodes*3-3+3;
% free_nodesDOF = [free_nodesDOF1, free_nodesDOF2, free_nodesDOF3];
% free_nodesDOF = sort(free_nodesDOF);
% free_u = openKKfoot(free_nodesDOF, free_nodesDOF)\...
%     (openQ(free_nodesDOF)-openKKfoot(free_nodesDOF,topNodesDOF3)*topMove*ones(length(topNodes),1));
% combined_u = zeros(length(openKKfoot),1);
% combined_u(free_nodesDOF) = free_u;
% combined_u(topNodesDOF3) = topMove;
% combined_u_allNode = zeros(length(KKfoot),1);
% combined_u_allNode(newNodeDOF) = combined_u;
% PlotDefoMesh([x', y', z'], wholeElem2n,10, [combined_u_allNode(1:3:end),combined_u_allNode(2:3:end),combined_u_allNode(3:3:end)])
% title('Trival test deformed mesh'); xlabel('x');ylabel('y');zlabel('z');
%----------------- trival test end--------------------------------------
%----------------- trival test 2--------------------------------------
% leftNodes = find(ox == -20);
% leftNodesDOF1 = leftNodes*3-3+1;
% leftNodesDOF2 = leftNodes*3-3+2;
% leftNodesDOF3 = leftNodes*3-3+3;
% leftNodesDOF = [leftNodesDOF1, leftNodesDOF2, leftNodesDOF3];
% leftNodesDOF = sort(leftNodesDOF);
% rightNodes = find(ox == 20);
% rightNodesDOF1 = rightNodes*3-3+1;
% rightNodesDOF2 = rightNodes*3-3+2;
% rightNodesDOF3 = rightNodes*3-3+3;
% rightNodesDOF = [rightNodesDOF1, rightNodesDOF2, rightNodesDOF3];
% rightNodesDOF = sort(rightNodesDOF);
% free_nodes = find(ox ~= -20&ox ~= 20);
% free_nodesDOF1 = free_nodes*3-3+1;
% free_nodesDOF2 = free_nodes*3-3+2;
% free_nodesDOF3 = free_nodes*3-3+3;
% free_nodesDOF = [free_nodesDOF1, free_nodesDOF2, free_nodesDOF3];
% free_nodesDOF = sort(free_nodesDOF);
% free_u = openKKfoot(free_nodesDOF, free_nodesDOF)\...
%     openQ(free_nodesDOF)*10;
% combined_u = zeros(length(openKKfoot),1);
% combined_u(free_nodesDOF) = free_u;
% combined_u_allNode = zeros(length(KKfoot),1);
% combined_u_allNode(newNodeDOF) = combined_u;
% PlotDefoMesh([x', y', z'], wholeElem2n,10, [combined_u_allNode(1:3:end),combined_u_allNode(2:3:end),combined_u_allNode(3:3:end)])
% title('Trival test deformed mesh'); xlabel('x');ylabel('y');zlabel('z');
%----------------- trival test 2 end--------------------------------------


ucat = zeros(length(openKKfoot),1);
ucat(groundNodeDOF1) = Uffx(oz == min(oz));
ucat(groundNodeDOF2) = Uffy(oz == min(oz));
ucat(groundNodeDOF3) = Uffz(oz == min(oz));
nodesPerLayer = length(bottom_nodes);
%% FLEXIBILITY MATRIX
no_connect = zeros(size(interNodesXYZ,1),1);
connect_xy = zeros(length(no_connect),8);
patchArea = zeros(size(interNodesXYZ,1),1);
elemAArray = zeros(length(elem2nInter),1);
for i=1:length(elem2nInter)
    meanX = mean(interNodesXYZ(elem2nInter(i,1:4),1));
    meanY = mean(interNodesXYZ(elem2nInter(i,1:4),2));
    vecA = interNodesXYZ(elem2nInter(i,3),:) - interNodesXYZ(elem2nInter(i,1),:);
    vecB = interNodesXYZ(elem2nInter(i,4),:) - interNodesXYZ(elem2nInter(i,2),:);
    elemA = norm(cross(vecA, vecB))/2;
    elemAArray(i) = elemA;
    for j=1:4
        node_no_old  = elem2nInter(i,j);
        no_connect(node_no_old) = no_connect(node_no_old)+1;
        connect_xy(node_no_old, no_connect(node_no_old)*2-1) = meanX;
        connect_xy(node_no_old, no_connect(node_no_old)*2-0) = meanY;
        patchArea(node_no_old) = patchArea(node_no_old) + elemA/4;
    end
end

dummyNodes = find(no_connect < 4);
patchAreaReal = patchArea;
patchAreaReal(dummyNodes) = [];
dummyElemBool = false(length(elem2nInter),1);
for i=1:length(elem2nInter)
    for j=1:4
        if sum(elem2nInter(i,j)==dummyNodes)~=0
            dummyElemBool(i) = true;
        end
    end
end

realNodes = unique(reshape(elem2nInter(~dummyElemBool,:), length(elem2nInter(~dummyElemBool,:))*4,1));
nus = nus*ones(size(realNodes));
Es = Es*ones(size(realNodes));
Gs = Es./2./(1+nus);
FLEX = computeFLEX16int_v4_test(interNodesXYZ,...
    elem2nInter(~dummyElemBool,:), connect_xy, Gs, nus);
FLEX = (FLEX + FLEX')*0.5;
FLEX = FLEX*1000000;
Ks_groundDOF = inv(FLEX);
ucat = ucat*1000;
openQ = openQ/1000;
openKKfoot = openKKfoot/1000000;
%----------------------trival 3
% u = (openKKfoot + Ks)\(openQ);
% u=u/1000;
% figure
% plot3(ox(oz==0), oy(oz==0), u(groundNodeDOF3),'x')
% xlabel('x')
% ylabel('y')
% ylabel('z')
% u_s = zeros(length(x),1);
% u_s(newNodeDOF) = u;
% PlotDefoMesh([x', y', z'], wholeElem2n,10, [u_s(1:3:end),u_s(2:3:end),u_s(3:3:end)])
% title('Deformed Mesh'); xlabel('x');ylabel('y');zlabel('z');
%----------------------
%----------------------elastic processing start------------------------------
% u = (openKKfoot + Ks)\(openQ+Ks*ucat);
% % u = (openKKfoot + Ks)\(openQ);
% u=u/1000;
% figure
% plot3(ox(oz==0), oy(oz==0), Uffz(oz==0,1),'x')
% hold on
% plot3(ox(oz==0), oy(oz==0), u(groundNodeDOF3),'x')
% xlabel('x')
% ylabel('y')
% ylabel('z')
% % frontNodes = find(oz==0&oy==-5);
% % figure
% % plot(ox(frontNodes), Uffz(frontNodes),'x')
% % hold on
% % plot(ox(frontNodes), u(frontNodes*3),'x')
% % title('Front facade vertical')
% % u_el_final = u;
% % rareNodes = find(oz==0&oy==5);
% % figure
% % plot(ox(rareNodes), Uffz(rareNodes),'x')
% % hold on
% % plot(ox(rareNodes), u(rareNodes*3,1),'x')
% % title('Rare facade vertical')
% % endNodes = find(oz==0 & ox ==10);
% % figure
% % plot(oy(endNodes), Uffz(endNodes),'x')
% % hold on
% % plot(oy(endNodes), u(endNodes),'x')
% % title('end facade vertical')
% eps = zeros(size(wholeElem2n, 1), 6);
% eps_principal = zeros(size(wholeElem2n, 1), 3);
% centroidXYZ = zeros(size(wholeElem2n, 1), 3);
% u_s = zeros(length(x),1);
% u_s(newNodeDOF) = u;
% % compute element strain
% for i = 1:size(wholeElem2n, 1)
%     xyzLocal = zeros(8, 3);
%     nGlobal = [];
% %     nGlobalVertical = [];
%     for j = 1:8
%         xyzLocal(j,:) = [x(wholeElem2n(i,j)), y(wholeElem2n(i,j)), z(wholeElem2n(i,j))];
%         nGlobal = [nGlobal,3*wholeElem2n(i,j)-2,3*wholeElem2n(i,j)-1,3*wholeElem2n(i,j)];
% %         nGlobalVertical = [nGlobalVertical,3*elem2n(i,j)];
%     end
%     [N, dNdx] = shape3d ([0 0 0],xyzLocal);
%     B = zeros(6,24);
%     B([1 4 6],1:3:24) = dNdx;
%     B([4 2 5],2:3:24) = dNdx;
%     B([6 5 3],3:3:24) = dNdx;
%     eps(i,:) = (B*u_s(nGlobal))';%xx,yy,zz,xy,yz,xz
%     eps_matrix = [eps(i,1),eps(i,4),eps(i,6);
%                   eps(i,4),eps(i,2),eps(i,5);
%                   eps(i,6),eps(i,5),eps(i,3)];
%     eps_principal(i,:) = eig(eps_matrix); 
%     centroidXYZ(i,:) = [N*xyzLocal(:,1),N*xyzLocal(:,2),N*xyzLocal(:,3)];
% end
% max_eps = max(eps);
% min_eps = min(eps);
% max_eps_principal = max(eps_principal);
% min_eps_principal = min(eps_principal);
% max_princpal_eps = [min_eps_principal(1),max_eps_principal(3)];
% PlotDefoMesh([x', y', z'], wholeElem2n,10, [u_s(1:3:end),u_s(2:3:end),u_s(3:3:end)])
% title('Deformed Mesh'); xlabel('x');ylabel('y');zlabel('z');
% hold on
% groundxN = 50;
% groundyN = 25;
% groundNodesX = linspace(-30, 30, groundxN);
% groundNodesY = linspace(-10, 10, groundyN);
% [gX,gY] = meshgrid(groundNodesX,groundNodesY);
% gZ = zeros(size(gX));
% for i=1:size(gX,1)
%     for j=1:size(gX,2)
% [gUffx(i,j),gUffy(i,j),gUffz(i,j)] = u_3D_comas_v2(Vltp,...
%         width_Kx, width_Ky, gX(i,j), gY(i,j), gZ(i,j), tunnel_delta, tunnel_z0,...
%         tunnel_ys, tunnel_yf, tunnel_d, 1);
%     end
% end
% surf(gX+gUffx,gY+gUffy,gZ+gUffz,'LineStyle','none');
% 
% PlotFieldonMesh([x', y', z'], wholeElem2n,u_s(3:3:end))
% title('u_z'); xlabel('x');ylabel('y');zlabel('z');
% PlotElemFieldonMesh([x', y', z'], wholeElem2n,eps_principal(:,1))
% title('Eps\_principal_1'); xlabel('x');ylabel('y');zlabel('z');
% PlotElemFieldonMesh([x', y', z'], wholeElem2n,eps_principal(:,3))
% title('Eps\_principal_3'); xlabel('x');ylabel('y');zlabel('z');
% PlotElemFieldonMesh([x', y', z'], wholeElem2n([1:1080, 2777:length(wholeElem2n)],:),eps_principal([1:1080, 2777:length(wholeElem2n)],3))
% % PlotElemFieldonMesh([x', y', z'], wholeElem2n([1:1272, 2969:length(wholeElem2n)],:),eps_principal([1:1272, 2969:length(wholeElem2n)],3))
% title('Eps\_principal_3 without foundation'); xlabel('x');ylabel('y');zlabel('z');
% PlotElemFieldonMesh([x', y', z'], wholeElem2n(rareElem,:),eps_principal(rareElem,3))
% title('Eps\_principal_3 rare wall'); xlabel('x');ylabel('y');zlabel('z');
% son&Cording Strain
% u=u_s;
% xa = [-24, -22.7, -21.7, -20.4, -19.6, -18.3, -17.3, -16,...
%       -14.7, -13.7, -12.4, -11.6, -10.3, -9.3, -8,...
%       -6.7, -5.7, -4.4, -3.6, -2.3, -1.3]+24;
% xa = [xa, xa+24];
% xb = [-22.7, -21.7, -20.4, -19.6, -18.3, -17.3, -16,...
%       -14.7, -13.7, -12.4, -11.6, -10.3, -9.3, -8,...
%       -6.7, -5.7, -4.4, -3.6, -2.3, -1.3, 0]+24;
% xb = [xb, xb+24];
% xc = xb;
% xd = xa;
% ya = ones(size(xa))*(-16);
% yb = ones(size(xa))*(-16);
% yc = ones(size(xa))*(-16);
% yd = ones(size(xa))*(-16);
% za = ones(size(xa))*(0.25);
% za(4) = 2.25;
% zb = ones(size(xa))*(0.25);
% zb(4) = 3.5;
% zc = ones(size(xa))*(3.25);
% zd = ones(size(xa))*(3.25);
% xa = [xa, xa];
% xb = [xb, xb];
% xc = [xc, xc];
% xd = [xd, xd];
% ya = [ya, ya+16];
% yb = [yb, yb+16];
% yc = [yc, yc+16];
% yd = [yd, yd+16];
% za = [za, za];
% zb = [zb, zb];
% zc = [zc, zc];
% zd = [zd, zd];
% theta = -26/180*pi;
% R = [cos(theta) -sin(theta); sin(theta), cos(theta)];
% for i=1:length(ox)
%     temp =  R*[ox(i); oy(i)];
%     xorth(i) = temp(1); yorth(i) = temp(2);
% end
% for i = 1:length(xa)
%     nodeInd = find(abs(xorth- xa(i))<1e-5 & abs(yorth- ya(i))<1e-5 & abs(oz-za(i))<1e-5);
%     ua = [u(nodeInd*3-2), u(nodeInd*3-1), u(nodeInd*3-0)];
%     nodeInd = find(abs(xorth- xb(i))<1e-5 & abs(yorth- yb(i))<1e-5 & abs(oz-zb(i))<1e-5);
%     ub = [u(nodeInd*3-2), u(nodeInd*3-1), u(nodeInd*3-0)];
%     nodeInd = find(abs(xorth- xc(i))<1e-5 & abs(yorth- yc(i))<1e-5 & abs(oz-zc(i))<1e-5);
%     uc = [u(nodeInd*3-2), u(nodeInd*3-1), u(nodeInd*3-0)];
%     nodeInd = find(abs(xorth- xd(i))<1e-5 & abs(yorth- yd(i))<1e-5 & abs(oz-zd(i))<1e-5);
%     ud = [u(nodeInd*3-2), u(nodeInd*3-1), u(nodeInd*3-0)];
%     openWallStrain(i,:) = sonAndCording3Dstrain(ua, ub, uc, ud, xb(i)-xa(i), 8, 1, 1);
% end
%----------------------elastic processing end------------------------------
%% --------------elastoplastic processing start------------------------------
freeNode = find(oz~=min(oz));
freeNodeDOF = freeNode*3-3+1;
freeNodeDOF = [freeNodeDOF, freeNode*3-3+2];
freeNodeDOF = [freeNodeDOF, freeNode*3-3+3];
freeNodeDOF = sort(freeNodeDOF);
K11 = openKKfoot(groundNodeDOF,groundNodeDOF);
K12 = openKKfoot(groundNodeDOF,freeNodeDOF);
K21 = openKKfoot(freeNodeDOF,groundNodeDOF);
K22 = openKKfoot(freeNodeDOF,freeNodeDOF);
K22invK21 = K22\K21;
Kcon = K11-K12*(K22\K21);
Qcon = K12*(K22\openQ(freeNodeDOF));
Kstar = inv(diag(diag(FLEX)));
Lstar = FLEX - diag(diag(FLEX));
LL = FLEX;
%%
Stiffness = Kcon + Kstar;
n_iter = 2000;
uinc = zeros(length(Qcon),1);
uip = zeros(length(Qcon),1);
ucap = zeros(length(Qcon),1);
uinc_prev = zeros(length(Qcon),1);
uip_prev = zeros(length(Qcon),1);
ucap_prev = zeros(length(Qcon),1);
duip = zeros(length(Qcon),1);
ducap = zeros(length(Qcon),1);
ducat = zeros(length(Qcon),1);
P_el = (-Qcon+openQ(groundNodeDOF)); 
dP = (-Qcon+openQ(groundNodeDOF))/num_incr_P;
beta = 4;
nodesPerLayer = length(groundNodeDOF)/3;
% for i = 1:num_incr_P
% iter = 0;
% err = 0;
% residual = 0;
% while iter < n_iter
%     iter = iter + 1;
% %     if iter == 2000
% %         "iter = 2000"
% %     end
%     du = Stiffness\(dP+Kstar*ducap+Kstar*ducat+Kstar*duip);
%     f_prev = (i-1)*(-Qcon+openQ(groundNodeDOF))/num_incr_P - Kcon*uinc_prev;
%     f_prev_v = f_prev(3:3:3*nodesPerLayer);
%     f_prev_x = f_prev(1:3:3*nodesPerLayer);
%     f_prev_y = f_prev(2:3:3*nodesPerLayer);
%     uinc = uinc_prev + du;
%     f_curr = i*(-Qcon+openQ(groundNodeDOF))/num_incr_P - Kcon*uinc;
%     f_curr_v = f_curr(3:3:3*nodesPerLayer);
%     df = dP-Kcon*du;
%     df_v = df(3:3:3*nodesPerLayer);
%     df_v(f_curr_v>0) = 0 - f_prev_v(f_curr_v>0);
%     f_v_capacity = -1*(100+8*(23-8)*0/2)*2/(1-0.425) * patchAreaReal;
%     df_v(f_curr_v<f_v_capacity) = f_v_capacity(f_curr_v<f_v_capacity) - f_prev_v(f_curr_v<f_v_capacity);
%     f_curr_v = df_v + f_prev_v;
%     f_h_max = max(-f_curr_v * mu,0);
%     f_h_x = f_curr(1:3:3*nodesPerLayer);
%     f_h_y = f_curr(2:3:3*nodesPerLayer);
%     f_h_curr = sqrt(f_h_x.^2 + f_h_y.^2);
%     slide_dof = abs(f_h_curr)>abs(f_h_max);
%     f_h_x_before = f_h_x;
%     f_h_y_before = f_h_y;
%     f_h_x(slide_dof) = f_h_x(slide_dof).*abs(f_h_max(slide_dof)./f_h_curr(slide_dof));
%     f_h_y(slide_dof) = f_h_y(slide_dof).*abs(f_h_max(slide_dof)./f_h_curr(slide_dof));
%     df_h_x = df(1:3:3*nodesPerLayer);
%     df_h_y = df(2:3:3*nodesPerLayer);
%     df_h_x(slide_dof)=f_h_x(slide_dof)-f_prev_x(slide_dof);
%     df_h_y(slide_dof)=f_h_y(slide_dof)-f_prev_y(slide_dof);
%     df(1:3:3*nodesPerLayer) = df_h_x;
%     df(2:3:3*nodesPerLayer) = df_h_y;
%     df(3:3:3*nodesPerLayer) = df_v;
% %     CI = (Lstar*(dP-KKfoot*du) + ducat + duip)./du;
%     du(abs(du)<10e-8) = 0;%% attantion
%     du_prim = (LL*df + ducat + duip);
%     du_prim(abs(du_prim)<10e-8) = 0;
%     CI = du_prim./du;
%     CI = CI(1:3*nodesPerLayer);
%     CI(isnan(CI(:,1)),:)=[1];
%     duip_old = duip;
%     duip = du-(LL-Lstar)*df-ducap-ducat;
%     ducap = (Lstar*df+beta*ducap)/(1+beta);
%     
% %     duip = (duip+duip_old)/2;
%     duip_prim = (duip-duip_old);
%     duip_prim(abs(du_prim)<10e-8) = 0;
%     perr = duip_prim./duip;
%     perr = perr(1:3*nodesPerLayer);
%     perr(isnan(perr(:,1)),:)=[0];
%     perr = norm(duip-duip_old);
%     err(iter,1) = max(abs(CI-1));
%     err(iter,2) = max(abs(perr));
%     uip = uip_prev + duip;
%     ucap = ucap_prev + ducap;
% %     residual(iter,4) = norm((Kcon + Kstar + Kstar * Lstar * Kcon)*uinc-...
% %         (Kstar*uip + Kstar * ducat*i + dP*i  +Kstar * Lstar * dP*i));
% %     u_bottom = uinc;
% %     u_solid = K22\(i*openQ(freeNodeDOF)/num_incr_P-K21*u_bottom);
% %     u_total = [u_bottom;u_solid];
% %     u_s_el = u_total-ucat*i/num_incr_TUNNEL*0;
% %     u_s_el(1:3*nodesPerLayer) = u_s_el(1:3*nodesPerLayer)-uip;
% %     residual(iter,2) = norm(KKfoot * u_total + Ks * u_s_el - ...
% %         i * Q/num_incr_P);
%     residual(iter,1) = norm(Stiffness*uinc - (dP*i+Kstar*ucap+Kstar*ducat+Kstar*uip));
% %     react_ep = P_el - Kcon*uinc;
%     react_ep = dP*i -Kcon*uinc;
%     react_ep_3dof = react_ep(3:3:3*nodesPerLayer);
%     residual(iter,3)=sum(react_ep_3dof>1);
%     if max(abs(CI-1))<0.05 && max(abs(perr))<0.05 && residual(iter,1)<0.00001
%         uinc_prev = uinc; 
%         uip_prev = uip;
%         ucap_prev = ucap;
%         break
%     end
%     if residual(iter,1)<10^(-6)
%        uinc_prev = uinc; 
%         uip_prev = uip;
%         ucap_prev = ucap;
%         break
%     end
% %     if iter>1 && residual(iter,1)>residual(iter-1,1) && residual(iter,1)<100
% %        uinc_prev = uinc; 
% %         uip_prev = uip;
% %         ucap_prev = ucap;
% %         break
% %     end
%     if iter ==n_iter
%         uinc_prev = uinc;
%         uip_prev = uip;
%         ucap_prev = ucap;
%     end
% end
% end
% Reaction1 = P_el - Kcon*uinc;
% Reaction1_3DOF = Reaction1(3:3:end);
% Support_check1 = sum(Reaction1_3DOF>1);
% Residual1 = (Kcon + Kstar)*uinc - (P_el + Kstar*ucap + ...
%     Kstar*uip);
% Residual1_2 = (Kcon + Kstar)*uinc - (P_el + Kstar*(Lstar*Reaction1) + ...
%     Kstar*uip);
% norm(Residual1)
%%
%-------------use elastic self weigth solution
% tic
% % u_P_el = (Kcon+Ks(groundNodeDOF, groundNodeDOF))\P_el;
u_P_el = (Kcon+Ks_groundDOF)\P_el;
% % u_P_eltest = (Kcon+Ks_groundDOFtest)\P_el;
% toc
figure
hold on 
plot3(ox(oz==min(z)), oy(oz==min(z)), Uffz(oz==min(z)),'x')
plot3(ox(oz==min(z)), oy(oz==min(z)), u_P_el(3:3:end)/1000,'x')
% plot3(ox(oz==min(z)), oy(oz==min(z)), u_P_eltest(3:3:end)/1000,'x')
xlabel('x'); ylabel('y'); zlabel('z (m)')
view(3)
uinc = u_P_el;
uip = zeros(size(uinc));
duip = zeros(size(uinc));
ducap = zeros(size(uinc));
uip_prev = zeros(size(uinc));
ucap_prev = Lstar*(P_el-Kcon*uinc);
%--------------------------------------------------------------------------

u_P = uinc;
uinc_prev = uinc;
dP = zeros(size(Qcon));
ducat = ucat(groundNodeDOF)/num_incr_TUNNEL;
% n_iter = 2000;
for i = 1:num_incr_TUNNEL
iter = 0;
err = 0;
residual = 0;
duip = zeros(size(uinc));
ducap = zeros(size(uinc));
while iter < n_iter
    iter = iter + 1;
%     tic
    du = Stiffness\(dP+Kstar*ducap+Kstar*ducat+Kstar*duip);
%     toc
    f_prev = (-Qcon+openQ(groundNodeDOF)) - Kcon*uinc_prev;
    f_prev_v = f_prev(3:3:3*nodesPerLayer);
    f_prev_x = f_prev(1:3:3*nodesPerLayer);
    f_prev_y = f_prev(2:3:3*nodesPerLayer);
    uinc = uinc_prev + du;
    f_curr = (-Qcon+openQ(groundNodeDOF)) - Kcon*uinc;
    f_curr_v = f_curr(3:3:3*nodesPerLayer);
    df = dP-Kcon*du;
    df_v = df(3:3:3*nodesPerLayer);
    df_v(f_curr_v>0) = 0 - f_prev_v(f_curr_v>0);
    f_v_capacity = -1*(100+8*(23-8)*0/2)*2/(1-0.425) * patchAreaReal;
%     f_v_capacity = -inf;
%     f_v_capacity = f_v_capacity';
    df_v(f_curr_v<f_v_capacity) = f_v_capacity(f_curr_v<f_v_capacity) - f_prev_v(f_curr_v<f_v_capacity);
    f_curr_v = df_v + f_prev_v;
    f_h_max = max(-f_curr_v * mu,0);
    f_h_x = f_curr(1:3:3*nodesPerLayer);
    f_h_y = f_curr(2:3:3*nodesPerLayer);
    f_h_curr = sqrt(f_h_x.^2 + f_h_y.^2);
    slide_dof = abs(f_h_curr)>abs(f_h_max);
    f_h_x(slide_dof) = f_h_x(slide_dof).*abs(f_h_max(slide_dof)./f_h_curr(slide_dof));
    f_h_y(slide_dof) = f_h_y(slide_dof).*abs(f_h_max(slide_dof)./f_h_curr(slide_dof));
    df_h_x = df(1:3:3*nodesPerLayer);
    df_h_y = df(2:3:3*nodesPerLayer);
    df_h_x(slide_dof)=f_h_x(slide_dof)-f_prev_x(slide_dof);
    df_h_y(slide_dof)=f_h_y(slide_dof)-f_prev_y(slide_dof);
    df(1:3:3*nodesPerLayer) = df_h_x;
    df(2:3:3*nodesPerLayer) = df_h_y;
    df(3:3:3*nodesPerLayer) = df_v;
%     CI = (Lstar*(dP-KKfoot*du) + ducat + duip)./du;
    du(abs(du)<10e-8) = 0;%% attention
    du_prim = (LL*df + ducat + duip);
    du_prim(abs(du_prim)<10e-8) = 0;
    CI = du_prim./du;
    CI = CI(1:3*nodesPerLayer);
    CI(isnan(CI(:,1)),:)=[1];
    duip_old = duip;
    
    duip = du-(LL-Lstar)*df-ducap-ducat;
    ducap = (Lstar*df+beta*ducap)/(1+beta);
    duip(abs(duip)<10e-8) = 0;%% attantion
    duip_prim = (duip-duip_old);
    duip_prim(abs(du_prim)<10e-8) = 0;
    perr = duip_prim./duip;
    perr = perr(1:3*nodesPerLayer);
    perr(isnan(perr(:,1)),:)=[0];
    err(iter,1) = max(abs(CI-1));
    err(iter,2) = max(abs(perr));
    uip = uip_prev + duip;
    ucap = ucap_prev + ducap;
%     residual(iter,4) = norm(Stiffness*uinc-(Kstar*uip +...
%         Kstar * ducat*i + P_el + ...
%         Kstar * Lstar * (P_el-Kcon*uinc)));
%     u_bottom = uinc;
%     u_solid = K22\(Q(freeNodeDOF)-K21*u_bottom);
%     u_total = [u_bottom;u_solid];
%     u_s_el = u_total-ucat*i/num_incr_TUNNEL;
%     u_s_el(1:3*nodesPerLayer) = u_s_el(1:3*nodesPerLayer)-uip;
%     residual(iter,2) = norm(KKfoot * u_total + Ks * u_s_el - Q);
    react_ep = P_el-Kcon*uinc;
    react_ep_3dof = react_ep(3:3:3*nodesPerLayer);
    residual(iter,3)=sum(react_ep_3dof>1);
    residual(iter,1) = norm(Stiffness*uinc - (P_el+Kstar*ucap+Kstar*ducat*i+...
        Kstar*uip));
    if iter == 1
        delta_residual = residual(iter,1)-0;
    else
        delta_residual = residual(iter,1)-residual(iter-1,1);
    end
%     if max(abs(CI-1))<0.05 && max(abs(perr))<0.05 && delta_residual<0.01 &&...
%             residual(iter,1)<1
    if max(abs(CI-1))<0.05 && max(abs(perr))<0.05 && residual(iter,1)<10^(-5)
        uinc_prev = uinc; 
        uip_prev = uip;
        ucap_prev = ucap;
        break
    end
    if residual(iter,1)<10^(-7)
       uinc_prev = uinc; 
        uip_prev = uip;
        ucap_prev = ucap;
        break
    end 
%     if iter>1 && residual(iter,1)>residual(iter-1,1) 
%        uinc_prev = uinc; 
%         uip_prev = uip;
%         ucap_prev = ucap;
% %         break
%     end
    if iter ==n_iter
        uinc_prev = uinc;
        uip_prev = uip;
        ucap_prev = ucap;
    end
end
end
toc
%% plot
u=uinc;
figure
plot3(ox(oz==min(oz)), oy(oz==min(oz)), Uffz(oz==min(oz),1),'x')
hold on
plot3(ox(oz==min(oz)), oy(oz==min(oz)), u(3:3:end)/1000,'s')
plot3(ox(oz==min(oz)), oy(oz==min(oz)), u_P(3:3:end)/1000,'^')
legend('Greenfield','Building final','Building selfweight')
xlabel('x')
ylabel('y')
ylabel('z')

% Final state:
u_bottom_final = u;
u_solid_final = K22\(openQ(freeNodeDOF)-K21*u_bottom_final);
u_total_final = zeros(length(ox),1);
u_total_final(groundNodeDOF) = u_bottom_final;
u_total_final(freeNodeDOF) = u_solid_final;
u_total_final = u_total_final/1000;
u_full_node_final = zeros(length(x),1);
u_full_node_final(newNodeDOF) = u_total_final;
frontNodes = find(oz==0&oy==-5);
% figure
% plot(ox(frontNodes), Uffz(frontNodes),'x')
% hold on
% plot(ox(frontNodes), u_total_final(frontNodes*3),'x')
% title('Front facade vertical')
% PlotDefoMesh([x', y', z'], wholeElem2n,10, [u_full_node_final(1:3:end),u_full_node_final(2:3:end),u_full_node_final(3:3:end)])
% title('Deformed Mesh'); xlabel('x');ylabel('y');zlabel('z');
% PlotFieldonMesh([x', y', z'], wholeElem2n,u_full_node_final(3:3:end))
% title('u_z'); xlabel('x');ylabel('y');zlabel('z');
% for i = 1:size(wholeElem2n, 1)
%     xyzLocal = zeros(8, 3);
%     nGlobal = [];
% %     nGlobalVertical = [];
%     for j = 1:8
%         xyzLocal(j,:) = [x(wholeElem2n(i,j)), y(wholeElem2n(i,j)), z(wholeElem2n(i,j))];
%         nGlobal = [nGlobal,3*wholeElem2n(i,j)-2,3*wholeElem2n(i,j)-1,3*wholeElem2n(i,j)];
% %         nGlobalVertical = [nGlobalVertical,3*elem2n(i,j)];
%     end
%     [N, dNdx] = shape3d ([0 0 0],xyzLocal);
%     B = zeros(6,24);
%     B([1 4 6],1:3:24) = dNdx;
%     B([4 2 5],2:3:24) = dNdx;
%     B([6 5 3],3:3:24) = dNdx;
%     eps(i,:) = (B*u_full_node_final(nGlobal))';%xx,yy,zz,xy,yz,xz
%     eps_matrix = [eps(i,1),eps(i,4),eps(i,6);
%                   eps(i,4),eps(i,2),eps(i,5);
%                   eps(i,6),eps(i,5),eps(i,3)];
%     eps_principal(i,:) = eig(eps_matrix); 
%     centroidXYZ(i,:) = [N*xyzLocal(:,1),N*xyzLocal(:,2),N*xyzLocal(:,3)];
% end


u_bottom = u-u_P;
u_solid = K22\(openQ(freeNodeDOF)-K21*u_bottom);
u_total = zeros(length(ox),1);
u_total(groundNodeDOF) = u_bottom;
u_total(freeNodeDOF) = u_solid;
u_total = u_total/1000;
eps = zeros(size(wholeElem2n, 1), 6);
eps_principal = zeros(size(wholeElem2n, 1), 3);
centroidXYZ = zeros(size(wholeElem2n, 1), 3);
u_s = zeros(length(x),1);
u_s(newNodeDOF) = u_total;
%%
figure
plot3(ox(oz==min(oz)), oy(oz==min(oz)), Uffz(oz==min(oz),1),'x')
hold on
plot3(ox(oz==min(oz)), oy(oz==min(oz)), u_bottom(3:3:end)/1000,'x')
xlabel('x')
ylabel('y')
ylabel('z')
legend('Greenfield','Building')

%% compute element strain
for i = 1:size(wholeElem2n, 1)
    xyzLocal = zeros(8, 3);
    nGlobal = [];
%     nGlobalVertical = [];
    for j = 1:8
        xyzLocal(j,:) = [x(wholeElem2n(i,j)), y(wholeElem2n(i,j)), z(wholeElem2n(i,j))];
        nGlobal = [nGlobal,3*wholeElem2n(i,j)-2,3*wholeElem2n(i,j)-1,3*wholeElem2n(i,j)];
%         nGlobalVertical = [nGlobalVertical,3*elem2n(i,j)];
    end
    [N, dNdx] = shape3d ([0 0 0],xyzLocal);
    B = zeros(6,24);
    B([1 4 6],1:3:24) = dNdx;
    B([4 2 5],2:3:24) = dNdx;
    B([6 5 3],3:3:24) = dNdx;
    eps(i,:) = (B*u_s(nGlobal))';%xx,yy,zz,xy,yz,xz
    eps_matrix = [eps(i,1),eps(i,4),eps(i,6);
                  eps(i,4),eps(i,2),eps(i,5);
                  eps(i,6),eps(i,5),eps(i,3)];
    eps_principal(i,:) = eig(eps_matrix); 
    centroidXYZ(i,:) = [N*xyzLocal(:,1),N*xyzLocal(:,2),N*xyzLocal(:,3)];
end
max_eps = max(eps);
min_eps = min(eps);
max_eps_principal = max(eps_principal);
min_eps_principal = min(eps_principal);
max_princpal_eps = [min_eps_principal(1),max_eps_principal(3)];
% PlotElemFieldonMesh([x', y', z'], wholeElem2n,eps_principal(:,1))
% title('Eps\_principal_1'); xlabel('x');ylabel('y');zlabel('z');
% PlotElemFieldonMesh([x', y', z'], wholeElem2n,eps_principal(:,3))
% title('Eps\_principal_3'); xlabel('x');ylabel('y');zlabel('z');
% largeStrainElem = find(eps_principal(:,3)>=quantile(eps_principal(:,3),0.99));

% PlotElemFieldonMesh([x', y', z'], wholeElem2n([1:1080, 2777:length(wholeElem2n)],:),eps_principal([1:1080, 2777:length(wholeElem2n)],3))
% PlotElemFieldonMesh([x', y', z'], wholeElem2n(rareElem,:),eps_principal(rareElem,3))
% title('Eps\_principal_3 rare wall'); xlabel('x');ylabel('y');zlabel('z');
%% Plot strain on deformaed mesh
figure
wallElemInd = [1:timberElemIndex(1)-1,...
    timberElemIndex(2)+1:length(wholeElem2n)];
PlotFieldonDefoMesh([x', y', z'], wholeElem2n(wallElemInd,:),10, [u_s(1:3:end),u_s(2:3:end),u_s(3:3:end)],eps_principal(wallElemInd,3)*1000000)
daspect([1 1 1])
% disp({'dv at corner 1 to 4: ', u_bottom(18*3), u_bottom(33*3), u_bottom(49*3), u_bottom(64*3)})
% disp({'dv total at corner 1 to 4: ', u_bottom_final(18*3), u_bottom_final(33*3), u_bottom_final(49*3), u_bottom_final(64*3)})
