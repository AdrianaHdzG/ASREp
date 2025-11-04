% demo1 of MasonMesh tool
% This demo create a mesh for a 20x10m masonary facade sitting on a 1m
% thick foundation
% Revision history:
%   Jinyan Zhao, jinyan_zhao@berkeley.edu, Jan 2023
% Cite As
%   
% 
clear 
close all
%% Input
% Define the coordinate of facade corners
baseCornersXYZ = [-40, -5.5, 0;
                  40, -5.5, 0;
                  40,  5.5, 0;
                 -40,  5.5, 0];


% Define the connectivity of the corners, the define of base wall has 
% to be counter clock wise             
baseWall = [1 2;
            2 3;
            3 4;
            4 1];
% Thickness of each wall, each wall can have different thickness       
WallThick = 0.2*ones(4,1);
% Thickness of foundation below each wall, all foundations can have
% different width
foundationWidth = 1*ones(size(baseWall,1),1);
% Height of the walls
WallHeight = 8;
% Height of the foundation
foundationHeight = 0.5;
% No opening in the facade in this demo
opening = [];
%--------------------Meshing parameters----------------------------
WallThickNoElem = 2; % Number of elements perpendicular to the facade plain
WallExterElemApproSize = 1;% Approximate element size in the facade plain
heightElemSize = 1; % Approximage facade element size in the vertical direction
foundationOutNoElem = 2; % Number of foundation element outside of the facade 
foundationInNoElem = 2; % Number of foundation element inside of the facade
foundationHeighNoElem = 2; % Number of foundation element in the vertical direction


%theta = 90;              % rotation angle in degrees
%Rz = [cosd(theta) -sind(theta) 0;
%      sind(theta)  cosd(theta) 0;
%      0            0           1];

% --- Apply rotation to all corner points ---

%baseCornersXYZ = baseCornersXYZ * Rz.';   % (N×3) * (3×3)



%% Generate facade nodes (wallNodesXYZ) and mesh (elem2n)
% nodes
WallLength = zeros(size(WallThick));
[baseExterNodesXYZ, WallNoElem, wallSegNo, segElemNo]= meshOnCurveWithFixNodes(baseCornersXYZ, baseWall, ...
    WallExterElemApproSize, []);
baseNodesXYZ = zeros(sum(WallNoElem)*(WallThickNoElem+1),3);
baseNodesXYZ(1:sum(WallNoElem),:) = baseExterNodesXYZ;
outerCornersXYZ = baseCornersXYZ;
for i=1:WallThickNoElem
    % find corners
    newCornersXYZ = offsetCurveIn(outerCornersXYZ, baseWall, WallThick/WallThickNoElem);
    % find mesh on curve
    baseNodesXYZ(1+sum(WallNoElem)*i:sum(WallNoElem)*(i+1),:) = ...
        meshOnCurveWithFixSegment(newCornersXYZ, baseCornersXYZ, baseWall, ...
            segElemNo, opening);    
    outerCornersXYZ = newCornersXYZ;
end
baseInnerCornersXYZ = newCornersXYZ;
if size(opening,1)>0
    openingZ = union(unique(opening(:,3)),unique(opening(:,5)));
    openingZ = [0; openingZ; WallHeight];
    openingZ = unique(openingZ);
else
    openingZ = [0; WallHeight];
end
endHeight = 0;
heightZ = 0;
for i=2:length(openingZ)
   headHeight = openingZ(i);
   segElemNo = floor((headHeight-endHeight)/heightElemSize)+1;
   if floor((headHeight-endHeight)/heightElemSize) == (headHeight-endHeight)/heightElemSize
       segElemNo = floor((headHeight-endHeight)/heightElemSize);
   end
   
   for k = 1:segElemNo
      heightZ = [heightZ; endHeight + k*(headHeight-endHeight)/segElemNo]; 
   end
   endHeight = heightZ(end);
end

% Elements 
elemPerLayerWall = sum(WallNoElem)*WallThickNoElem;
WallHeightNoElem = length(heightZ)-1;
elem2n = zeros(elemPerLayerWall*WallHeightNoElem, 8);
for j = 1:WallThickNoElem
    for i=1:sum(WallNoElem)
       elem2n(i+(j-1)*sum(WallNoElem),1) = i + (j-1)*sum(WallNoElem);
       elem2n(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-1)*sum(WallNoElem);
       elem2n(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       elem2n(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       if i==sum(WallNoElem)
           elem2n(i+(j-1)*sum(WallNoElem),1) = i + (j-1)*sum(WallNoElem);
           elem2n(i+(j-1)*sum(WallNoElem),2) = 1 + (j-1)*sum(WallNoElem);
           elem2n(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem);
           elem2n(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       end
    end
end
elem2n(1:elemPerLayerWall,5:8) = ...
    elem2n(1:elemPerLayerWall,1:4) + size(baseNodesXYZ,1);
wallNodesXYZ = zeros(size(baseNodesXYZ,1)*(WallHeightNoElem+1),3);
wallNodesXYZ(1:size(baseNodesXYZ,1),:) = baseNodesXYZ;
for i = 1:WallHeightNoElem-1
   elem2n(i*elemPerLayerWall+1:(i+1)*elemPerLayerWall,:) = ...
       elem2n((i-1)*elemPerLayerWall+1:i*elemPerLayerWall,:) + ...
       size(baseNodesXYZ,1);
   wallNodesXYZ(i*size(baseNodesXYZ,1)+1:(i+1)*size(baseNodesXYZ,1),1:2) =...
       baseNodesXYZ(:,1:2);
   wallNodesXYZ(i*size(baseNodesXYZ,1)+1:(i+1)*size(baseNodesXYZ,1),3) = ...
       heightZ(i+1);
end
wallNodesXYZ(WallHeightNoElem*size(baseNodesXYZ,1)+1:(WallHeightNoElem+1)*size(baseNodesXYZ,1),1:2) =...
       baseNodesXYZ(:,1:2);
wallNodesXYZ(WallHeightNoElem*size(baseNodesXYZ,1)+1:(WallHeightNoElem+1)*size(baseNodesXYZ,1),3) = ...
       heightZ(end);

PlotMesh(wallNodesXYZ,elem2n, 0)

%% Add foundation nodes and element to get wholeNodesXYZ and wholeElem2n
% nodes
foundationNodesXYZ = zeros(sum(WallNoElem)*(WallThickNoElem+1+foundationOutNoElem+foundationInNoElem) * ...
    (foundationHeighNoElem+1) - size(baseNodesXYZ,1),3);
foundationNodesXYZ(1:sum(WallNoElem),:) = offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
    wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),WallThick/WallThickNoElem,...
    (foundationWidth-WallThick)/2);

for i = 1:foundationOutNoElem-1
    foundationNodesXYZ(1+i*sum(WallNoElem):(i+1)*sum(WallNoElem),:) = ...
        offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
        wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),WallThick/WallThickNoElem,...
        (foundationWidth-WallThick)/2-i*(foundationWidth-WallThick)/2/foundationOutNoElem);
end
for i = 1:foundationInNoElem
    foundationNodesXYZ(1+(i-1+foundationOutNoElem)*sum(WallNoElem):(i-1+1+foundationOutNoElem)*sum(WallNoElem),:) = ...
        offsetWallNodes(wallNodesXYZ(1+sum(WallNoElem)*WallThickNoElem:sum(WallNoElem)+sum(WallNoElem)*WallThickNoElem,:),...
        wallNodesXYZ(1+sum(WallNoElem)*(WallThickNoElem-1):sum(WallNoElem)*WallThickNoElem,:),WallThick/WallThickNoElem,...
        (foundationWidth-WallThick)/2-(i-1)*(foundationWidth-WallThick)/2/foundationInNoElem);
end

foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)+1:...
    sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)+sum(WallNoElem)*foundationOutNoElem,:) = ...
    foundationNodesXYZ(1:sum(WallNoElem)*(0+foundationOutNoElem),:);

foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
    +sum(WallNoElem)*foundationOutNoElem+1:...
    sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
    +sum(WallNoElem)*foundationOutNoElem + size(baseNodesXYZ,1),:) = ...
    baseNodesXYZ;

foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
    +sum(WallNoElem)*foundationOutNoElem+size(baseNodesXYZ,1)+1:...
    sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
    +sum(WallNoElem)*foundationOutNoElem + size(baseNodesXYZ,1)...
    +sum(WallNoElem)*foundationInNoElem,:) = ...
    foundationNodesXYZ(sum(WallNoElem)*foundationOutNoElem+1:...
    sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem),:);


foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)+1:...
    sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)*2 + size(baseNodesXYZ,1),3) = ...
    0 - foundationHeight/foundationHeighNoElem;

foundationLayerNoNodes = sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem) + size(baseNodesXYZ,1);
for i = 1:foundationHeighNoElem-1
    foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +i*foundationLayerNoNodes + 1:...
        sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +(i+1)*foundationLayerNoNodes,:) = ...
    foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +(i-1)*foundationLayerNoNodes + 1:...
        sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +(i)*foundationLayerNoNodes,:);

    foundationNodesXYZ(sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +i*foundationLayerNoNodes + 1:...
        sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
        +(i+1)*foundationLayerNoNodes,3) = 0 - foundationHeight/foundationHeighNoElem*(i+1);
end

wholeNodesXYZ = [wallNodesXYZ;foundationNodesXYZ];

% Elements
elem2nFound = zeros((WallThickNoElem+foundationOutNoElem+foundationInNoElem)*...
    sum(WallNoElem) * foundationHeighNoElem,8);
for j= 1:foundationOutNoElem
    for i = 1:sum(WallNoElem)
       elem2nFound(i+(j-1)*sum(WallNoElem),1) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),5) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),6) = i+1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),7) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       if j==foundationOutNoElem
           elem2nFound(i+(j-1)*sum(WallNoElem),3) = i+1;
           elem2nFound(i+(j-1)*sum(WallNoElem),4) = i;
           elem2nFound(i+(j-1)*sum(WallNoElem),7) = i+1 + size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
               + sum(WallNoElem)*foundationOutNoElem;
           elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
               + sum(WallNoElem)*foundationOutNoElem;
       end
       if i==sum(WallNoElem)
           elem2nFound(i+(j-1)*sum(WallNoElem),1) = i + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),2) = 1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),5) = i + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),6) = 1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),7) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           if j==foundationOutNoElem
               elem2nFound(i+(j-1)*sum(WallNoElem),3) = 1;
               elem2nFound(i+(j-1)*sum(WallNoElem),4) = i;
               elem2nFound(i+(j-1)*sum(WallNoElem),7) = 1 + size(wallNodesXYZ,1)...
                   + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
                   + sum(WallNoElem)*foundationOutNoElem;
               elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + size(wallNodesXYZ,1)...
                   + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
                   + sum(WallNoElem)*foundationOutNoElem;
           end
       end  
    end
end

for j= 1+foundationOutNoElem:WallThickNoElem+foundationOutNoElem
    for i = 1:sum(WallNoElem)
       elem2nFound(i+(j-1)*sum(WallNoElem),1) = i+(j-foundationOutNoElem-1)*sum(WallNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-foundationOutNoElem-1)*sum(WallNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-foundationOutNoElem-1)*sum(WallNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-foundationOutNoElem-1)*sum(WallNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),5) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),6) = i+1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),7) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       if i==sum(WallNoElem)
           elem2nFound(i+(j-1)*sum(WallNoElem),1) = i + (j-foundationOutNoElem-1)*sum(WallNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),2) = 1 + (j-foundationOutNoElem-1)*sum(WallNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-foundationOutNoElem-1)*sum(WallNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-foundationOutNoElem-1)*sum(WallNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),5) = i + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),6) = 1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),7) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       end  
    end
end

for j= 1+foundationOutNoElem+WallThickNoElem:...
        WallThickNoElem+foundationOutNoElem+foundationInNoElem
    for i = 1:sum(WallNoElem)
       elem2nFound(i+(j-1)*sum(WallNoElem),1) = i+(j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
       elem2nFound(i+(j-1)*sum(WallNoElem),5) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),6) = i+1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),7) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
           + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
       if j==1+foundationOutNoElem+WallThickNoElem
           elem2nFound(i+(j-1)*sum(WallNoElem),1) = i + sum(WallNoElem)*WallThickNoElem;
           elem2nFound(i+(j-1)*sum(WallNoElem),2) = i + 1 + sum(WallNoElem)*WallThickNoElem;
           elem2nFound(i+(j-1)*sum(WallNoElem),5) = i + size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
               + sum(WallNoElem)*(foundationOutNoElem + WallThickNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),6) = i + 1 + size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
               + sum(WallNoElem)*(foundationOutNoElem + WallThickNoElem);
       end
       if i==sum(WallNoElem)
           elem2nFound(i+(j-1)*sum(WallNoElem),1) = i+(j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),2) = 1 + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-WallThickNoElem-1-1)*sum(WallNoElem)+ size(wallNodesXYZ,1);
           elem2nFound(i+(j-1)*sum(WallNoElem),5) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),6) = 1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),7) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
               + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
           if j==1+foundationOutNoElem+WallThickNoElem
               elem2nFound(i+(j-1)*sum(WallNoElem),1) = i + sum(WallNoElem)*WallThickNoElem;
               elem2nFound(i+(j-1)*sum(WallNoElem),2) = 1 + sum(WallNoElem)*WallThickNoElem;
               elem2nFound(i+(j-1)*sum(WallNoElem),5) = i + size(wallNodesXYZ,1)...
                   + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
                   + sum(WallNoElem)*(foundationOutNoElem + WallThickNoElem);
               elem2nFound(i+(j-1)*sum(WallNoElem),6) = 1 + size(wallNodesXYZ,1)...
                   + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
                   + sum(WallNoElem)*(foundationOutNoElem + WallThickNoElem);
           end
       end  
    end
end

for i = 2:foundationHeighNoElem
    elem2nFound(sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*(i-1)+1:...
        sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*i, 1:4) = ...    
    elem2nFound(sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*(i-2)+1:...
        sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*(i-1), 5:8);
    
    elem2nFound(sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*(i-1)+1:...
        sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*i, 5:8) = ...
    elem2nFound(sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*(i-1)+1:...
        sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem)*i, 1:4)...
        + sum(WallNoElem)*(foundationOutNoElem+WallThickNoElem+foundationInNoElem+1);    
end

elem2nFoundReorder = elem2nFound;
elem2nFoundReorder(:,1:4) = elem2nFound(:,5:8);
elem2nFoundReorder(:,5:8) = elem2nFound(:,1:4);
wholeElem2n = [elem2n; elem2nFoundReorder];
PlotMesh(wholeNodesXYZ, wholeElem2n,0)




%% Save

timeString = datestr(datetime('now'),30);
timberElemIndex = [1 1];
title = strcat('demo_mat_',timeString,'.mat');
%save(title,'wholeNodesXYZ','wholeElem2n', 'timberElemIndex')
