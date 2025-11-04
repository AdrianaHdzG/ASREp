% demo4 of MasonMesh tool (A multi-story building)
% This demo create a mesh for the building Kansler gata 10 in Oslo Norway
%
% Revision history:
%   Jinyan Zhao, jinyan_zhao@berkeley.edu, Jan 2023
% Cite As
%   
% 
clear 
close all
%% Input 
% Define the coordinate of facade corners 
baseCornersXYZ = [ 64.98, 39.60, 0;
                   61.34, 21.02, 0;
                   75.55, 08.31, 0;
                   83.51, 17.23, 0;
                   74.76, 25.06, 0;
                   76.35, 32.77, 0;
                   86.30, 30.75, 0;
                   87.15, 35.83, 0];
baseCornersXYZ(:,1) = baseCornersXYZ(:,1) + 598700;
baseCornersXYZ(:,2) = baseCornersXYZ(:,2) + 6642100;
% Define the connectivity of the corners, the define of base wall has 
% to be counter clock wise              
baseWall = [1 2;
            2 3;
            3 4;
            4 5;
            5,6;
            6,7;
            7,8;
            8,1];
% Thickness of each wall, each wall can have different thickness 
WallThick = 0.35*ones(8,1);
% Thickness of foundation below each wall, all foundations can have
% different width
foundationWidth = 0.7*ones(size(baseWall,1),1);
% Height of the walls
WallHeight = 13.8;
% Height of the foundation
foundationHeight = 2;
%--------------------Meshing parameters------------------------------------
WallThickNoElem = 2; % Number of elements perpendicular to the facade plain
WallExterElemApproSize = 1;% Approximate element size in the facade plain
heightElemSize = 1; % Approximage facade element size in the vertical direction
foundationOutNoElem = 1; % Number of foundation element outside of the facade
foundationInNoElem = 1; % Number of foundation element inside of the facade
foundationHeighNoElem = 2; % Number of foundation element in the vertical direction
%--------------------- input of opening------------------------------------
% opening is defined as wall_no, x1/y1 z1 x2/y2 z2 
opening = [];
opening = [opening; [1 1.63 1.89 2.75 3.71]];
opening = [opening; [1 4.85 1.89 5.97 3.71]];
opening = [opening; [1 7.27 1.89 8.39 3.71]];
opening = [opening; [1 10.19 1.89 11.31 3.71]];
opening = [opening; [1 12.71 1.89 13.83 3.71]];
opening = [opening; [1 16.03 1.89 17.15 3.71]];
for i = 1:3
   temp = opening(1:6,:);
   temp(:,3) = temp(:,3)+3.2*i;
   temp(:,5) = temp(:,5)+3.2*i;
   opening = [opening;temp];
end
opening2 = [2, 2.16, 1.89, 3.28, 3.71];
opening2 = [opening2; [2 6.26 1.89 7.38 3.71]];
opening2 = [opening2; [2 9.58 1.89 10.70 3.71]];
opening2 = [opening2; [2 11.6 1.89 12.72 3.71]];
opening2 = [opening2; [2 14.52 1.89 15.64 3.71]];
opening2 = [opening2; [2 16.0 0 18.41 4.2]];
for i = 1:3
   temp = opening2(1:6,:);
   temp(6,3) = 1.89;temp(6,4) = 17.12;temp(6,5) = 3.71;
   temp(:,3) = temp(:,3)+3.2*i;
   temp(:,5) = temp(:,5)+3.2*i;
   opening2 = [opening2;temp];
end
opening = [opening; opening2];
opening4 = [4, 0.65, 0, 3.06, 4.2];
opening4 = [opening4; [4 3.56 1.89 4.68 3.71]];
opening4 = [opening4; [4 5.61 1.89 6.73 3.71]];
opening4 = [opening4; [4 8.5 1.89 9.62 3.71]];
opening4 = [opening4; [4 10.12 1.89 11.24 3.71]];
for i = 1:3
   temp = opening4(1:5,:);
   temp(1,3) = 1.89;temp(1,2) = 1.94;temp(1,5) = 3.71;
   temp(:,3) = temp(:,3)+3.2*i;
   temp(:,5) = temp(:,5)+3.2*i;
   opening4 = [opening4;temp];
end
opening = [opening; opening4];
opening5 = [5, 0.5, 1.89, 1.62, 3.71];
opening5 = [opening5; [5 2.94 1.89 4.06 3.71]];
opening5 = [opening5; [5 5.0 1.89 6.12 3.71]];
opening5 = [opening5; [5 6.42 1.89 7.37 3.71]];
for i = 1:3
   temp = opening5(1:4,:);
   temp(:,3) = temp(:,3)+3.2*i;
   temp(:,5) = temp(:,5)+3.2*i;
   opening5 = [opening5;temp];
end
opening = [opening; opening5];
opening6 = [6, 0.5, 1.89, 1.45, 3.71];
opening6 = [opening6; [6 2.86 1.89 3.98 3.71]];
opening6 = [opening6; [6 7.90 1.89 9.02 3.71]];
for i = 1:3
   temp = opening6(1:3,:);
   temp(:,3) = temp(:,3)+3.2*i;
   temp(:,5) = temp(:,5)+3.2*i;
   opening6 = [opening6;temp];
end
opening = [opening; opening6];
%% Generate facade nodes (wallNodesXYZ) and mesh (elem2n)
% nodes
WallLength = zeros(size(WallThick));
WallOpenNo = zeros(size(WallThick));
[baseExterNodesXYZ, WallNoElem, wallSegNo, segElemNo]= meshOnCurveWithFixNodes(baseCornersXYZ, baseWall, ...
    WallExterElemApproSize, opening);
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
%% Remove open element
% delete any brick elements that sit where a window/door opening should be 
% → turning the solid wall into a wall with holes.
windowElemNo = [];
for i = 1:size(opening,1)
    wall_no = opening(i, 1);
    if wall_no>1
       exterNoX = sum(WallNoElem(1:(wall_no-1)))+1;
       innerNoX = sum(WallNoElem(1:(wall_no-1)))+1;
    else
       exterNoX = 0+1;
       innerNoX = 0+1;
    end
    wc1xyz = baseCornersXYZ(baseWall(wall_no,1),:);
    wc2xyz = baseCornersXYZ(baseWall(wall_no,2),:);
    for j = exterNoX:sum(WallNoElem(1:wall_no))-1
        for k = 1:WallHeightNoElem
            for l = 1:WallThickNoElem
                elem_no = j + (k-1)*sum(WallNoElem)*WallThickNoElem + ...
                    (l-1)*sum(WallNoElem);
                meanX = mean(wallNodesXYZ(elem2n(elem_no,:),1));
                meanY = mean(wallNodesXYZ(elem2n(elem_no,:),2));
                meanZ = mean(wallNodesXYZ(elem2n(elem_no,:),3));
                projXYZ = projectOnLine([meanX,meanY,meanZ], wc1xyz, wc2xyz);
                projXYZ(3)=0;
                if opening(i,2)<dist(projXYZ, wc1xyz) && ...
                        opening(i,4)>dist(projXYZ, wc1xyz) && ...
                        opening(i,3)<meanZ && ...
                        opening(i,5)>meanZ
                    windowElemNo = [windowElemNo,elem_no];
                end
            end 
        end
    end       
end
elem2n(windowElemNo,:) = [];
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
PlotMesh(wholeNodesXYZ(:,1:2), wholeElem2n,0)
