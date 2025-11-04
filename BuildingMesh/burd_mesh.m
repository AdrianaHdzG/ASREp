clear 
close all
%% Input 
baseCornersXYZ = [-20, -5, 0;
                  20, -5, 0;
                  20,  5, 0;
                 -20,  5, 0];

% The define of base wall has to be counter clock wise 
% which corners are connected into walls 
baseWall = [1 2;
            2 3;
            3 4;
            4 1];

WallThick = 0.215*ones(4,1);
WallThickNoElem = 2;
WallExterElemApproSize = 1;
FootThick = 1*ones(4,1);
WallHeight = 8.5;
heightElemSize = 1;

% Input of opening
% opening is defined as wall_no, x1_or_y1 z1 x2_or_y2 z2 
% --- Example: add window on wall 1, 2.25–3.75 m along the wall, 1–3 m high (z)
opening = [];
opening = [opening; [1 2.25 1 3.75 3]];
opening = [opening; [1 2.25 5 3.75 7]];
opening = [opening; [1 5.25 1 6.75 3]];
opening = [opening; [1 5.25 5 6.75 7]];
opening = [opening; [1 9 0 11 3]];
opening = [opening; [1 9.25 5 10.75 7]];
opening = [opening; [1 13.25 1 14.75 3]];
opening = [opening; [1 13.25 5 14.75 7]];
opening = [opening; [1 16.25 1 17.75 3]];
opening = [opening; [1 16.25 5 17.75 7]];

% --- Mirror the same pattern to the right half of wall 1
opening2 = opening;
opening2(:,2) = opening2(:,2)+20;
opening2(:,4) = opening2(:,4)+20;
opening = [opening; opening2];

% --- Copy entire pattern from wall 1 to wall 3 (opposite side)
opening2 = opening;
opening2(:,1) = 3;
opening = [opening; opening2];

% --- Add a few extra small windows explicitly on wall 3
opening = [opening; [3 2.25 1 3.75 3]];
opening = [opening; [3 2.25 5 3.75 7]];
opening = [opening; [3 5.25 1 6.75 3]];
opening = [opening; [3 5.25 5 6.75 7]];


% --- Step 5: extra window near far-right end of wall 1
opening = [opening; [1 36.25 1 37.75 3]];

% --- Zero-height "split markers" at mid-wall (for mesh alignment)
% These are not real openings; they just force nodes to exist
% at the wall’s mid-point (20 m) and at ±half wall-thickness offsets
opening = [opening; [1 20-0.215/2 0 20+0.215/2 0]];
opening = [opening; [1 20 0 20 0]];
opening = [opening; [3 20-0.215/2 0 20+0.215/2 0]];
opening = [opening; [3 20 0 20 0]];

%--------------------Foundation----------------------------
foundationHeight = 0.5;
foundationWidth = 1*ones(size(baseWall,1),1);
foundationOutNoElem = 2;
foundationInNoElem = 2;
foundationHeighNoElem = 2;


%%--- CONFIG: rotation ---
applyRotation = true;   % set to true to rotate mesh at the end
theta = 90;         % rotation angle in degrees (about Z)
Rz = [cosd(theta) -sind(theta) 0;  % rotation matrix about Z
      sind(theta)  cosd(theta) 0;
      0            0           1];


%% Generate nodes
% WallExterElemTrueSize = zeros(size(WallThick));
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
%     figure    
%     plot(baseNodesXYZ(1+sum(WallNoElem)*i:sum(WallNoElem)*(i+1),1),baseNodesXYZ(1+sum(WallNoElem)*i:sum(WallNoElem)*(i+1),2),'x');
%     daspect([1 1 1])    
    outerCornersXYZ = newCornersXYZ;
end
baseInnerCornersXYZ = newCornersXYZ;

% plot(baseNodesXYZ(:,1),baseNodesXYZ(:,2),'x');
% daspect([1 1 1])
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

% wallNodesXYZ = zeros(size(baseNodesXYZ,1)*length(heightZ),3);
% for i = 1:length(heightZ)
%     wallNodesXYZ(1+(i-1)*size(baseNodesXYZ,1):i*size(baseNodesXYZ,1),:) = ...
%         baseNodesXYZ;
%     wallNodesXYZ(1+(i-1)*size(baseNodesXYZ,1):i*size(baseNodesXYZ,1),3) = ...
%         heightZ(i);
% end

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
PlotMesh(baseNodesXYZ(:,1:2),elem2n(1:elemPerLayerWall,1:4), 0)
% PlotMesh(wallNodesXYZ,elem2n, 0)


%% Remove open element
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
PlotMesh(wallNodesXYZ,elem2n, 0)

%% Add footing
startFoundationElem = length(elem2n)+1;
foundationNodesXYZ = zeros(sum(WallNoElem)*(WallThickNoElem+1+foundationOutNoElem+foundationInNoElem) * ...
    (foundationHeighNoElem+1) - size(baseNodesXYZ,1),3);
foundationNodesXYZ(1:sum(WallNoElem),:) = offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
    wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),WallThick/WallThickNoElem,...
    (foundationWidth-WallThick)/2);
for i = 1:foundationOutNoElem-1
%     foundationNodesXYZ(1+i*sum(WallNoElem):(i+1)*sum(WallNoElem),:) = ...
%         offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
%         wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),WallThick/WallThickNoElem,...
%         (foundationWidth-WallThick)/2-i*(foundationWidth-WallThick)/2/foundationOutNoElem, WallNoElem);
    foundationNodesXYZ(1+i*sum(WallNoElem):(i+1)*sum(WallNoElem),:) = ...
        offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
        wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),WallThick/WallThickNoElem,...
        (foundationWidth-WallThick)/2-i*(foundationWidth-WallThick)/2/foundationOutNoElem);
end
for i = 1:foundationInNoElem
%     foundationNodesXYZ(1+(i-1+foundationOutNoElem)*sum(WallNoElem):(i-1+1+foundationOutNoElem)*sum(WallNoElem),:) = ...
%         offsetWallNodes(wallNodesXYZ(1+sum(WallNoElem)*WallThickNoElem:sum(WallNoElem)+sum(WallNoElem)*WallThickNoElem,:),...
%         wallNodesXYZ(1+sum(WallNoElem)*(WallThickNoElem-1):sum(WallNoElem)*WallThickNoElem,:),WallThick/WallThickNoElem,...
%         (foundationWidth-WallThick)/2-(i-1)*(foundationWidth-WallThick)/2/foundationInNoElem, WallNoElem);
%     foundationNodesXYZ(1+(i-1+foundationOutNoElem)*sum(WallNoElem):(i-1+1+foundationOutNoElem)*sum(WallNoElem),:) = ...
%         offsetWallNodes(wallNodesXYZ(1+sum(WallNoElem)*WallThickNoElem:sum(WallNoElem)+sum(WallNoElem)*WallThickNoElem,:),...
%         wallNodesXYZ(1+sum(WallNoElem)*(WallThickNoElem-1):sum(WallNoElem)*WallThickNoElem,:),WallThick/WallThickNoElem,...
%         (foundationWidth-WallThick)/2-(i-1)*(foundationWidth-WallThick)/2/foundationInNoElem);
    foundationNodesXYZ(1+(i-1+foundationOutNoElem)*sum(WallNoElem):(i-1+1+foundationOutNoElem)*sum(WallNoElem),:) = ...
        offsetWallNodes(wallNodesXYZ(1+sum(WallNoElem)*WallThickNoElem:sum(WallNoElem)+sum(WallNoElem)*WallThickNoElem,:),...
        wallNodesXYZ(1+sum(WallNoElem)*(WallThickNoElem-1):sum(WallNoElem)*WallThickNoElem,:),WallThick/WallThickNoElem,...
        (i)*(foundationWidth-WallThick)/2/foundationInNoElem);
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
endFoundationElem = length(wholeElem2n);


%% Add internal wall
interWallCorners = [0 -5+0.215 0;
                    0  5-0.215 0];
interWall = [1 2];
interOpening = [];
InterNodesXYZ = [];
[InterCentralNodesXYZ, InterWallNoElem, wallSegNo, segElemNo]= meshOnCurveWithFixNodes(interWallCorners, interWall, ...
    WallExterElemApproSize, interOpening);
InterCentralNodesXYZ(end,:) = [];
InterNodesXYZ = [InterNodesXYZ; InterCentralNodesXYZ];


interWallCorners = [-0.215/2 -5+0.215 0;
                    -0.215/2  5-0.215 0];
[InterCentralNodesXYZ, InterWallNoElem, wallSegNo, segElemNo]= meshOnCurveWithFixNodes(interWallCorners, interWall, ...
    WallExterElemApproSize, interOpening);
InterCentralNodesXYZ(end,:) = [];
InterNodesXYZ = [InterNodesXYZ; InterCentralNodesXYZ];

interWallCorners = [+0.215/2 -5+0.215 0;
                    +0.215/2  5-0.215 0];
[InterCentralNodesXYZ, InterWallNoElem, wallSegNo, segElemNo]= meshOnCurveWithFixNodes(interWallCorners, interWall, ...
    WallExterElemApproSize, interOpening);
InterCentralNodesXYZ(end,:) = [];
InterNodesXYZ = [InterNodesXYZ; InterCentralNodesXYZ];
temp = InterNodesXYZ;



for i = 2:length(heightZ)
    temp(:,3) = heightZ(i);
    InterNodesXYZ = [InterNodesXYZ; temp];
end
    
wholeNodesXYZ = [wholeNodesXYZ; InterNodesXYZ];

interElem2n = zeros(InterWallNoElem*WallThickNoElem*WallHeightNoElem,8);


for i = 1:InterWallNoElem
    for j = 1:WallHeightNoElem
       if i==1
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,1) = ...
              find(abs(wholeNodesXYZ(:,1)-(-0.215/2))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j) &wholeNodesXYZ(:,2)== (-5+0.215/2));
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,2) = ...
              find(abs(wholeNodesXYZ(:,1)-       (0))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j) &wholeNodesXYZ(:,2)== (-5+0.215/2));
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,3) = ...
              i + (j-1)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,4) = ...
              i + (j-1)*3*(InterWallNoElem-1) + (InterWallNoElem-1)+ size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,5) = ...
              find(abs(wholeNodesXYZ(:,1)-(-0.215/2))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j+1) &wholeNodesXYZ(:,2)== (-5+0.215/2));
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,6) = ...
              find(abs(wholeNodesXYZ(:,1)-       (0))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j+1) &wholeNodesXYZ(:,2)== (-5+0.215/2));
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,7) = ...
              i + (j)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,8) = ...
              i + (j)*3*(InterWallNoElem-1) + (InterWallNoElem-1)+ size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,1) = ...
              find(abs(wholeNodesXYZ(:,1)-       (0))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j) &wholeNodesXYZ(:,2)== (-5+0.215/2));
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,2) = ...
              find(abs(wholeNodesXYZ(:,1)-(+0.215/2))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j) &wholeNodesXYZ(:,2)== (-5+0.215/2));
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,3) = ...
              i + (j-1)*3*(InterWallNoElem-1)+ (InterWallNoElem-1)*2 + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,4) = ...
              i + (j-1)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,5) = ...
              find(abs(wholeNodesXYZ(:,1)-       (0))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j+1) &wholeNodesXYZ(:,2)== (-5+0.215/2));
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,6) = ...
              find(abs(wholeNodesXYZ(:,1)-(+0.215/2))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j+1) &wholeNodesXYZ(:,2)== (-5+0.215/2));
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,7) = ...
              i + (j)*3*(InterWallNoElem-1) +(InterWallNoElem-1)*2+size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,8) = ...
              i + (j)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
       elseif i==InterWallNoElem
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,1) = ...
              i-1 + (j-1)*3*(InterWallNoElem-1) + (InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,2) = ...
              i-1 + (j-1)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,3) = ...
              find(abs(wholeNodesXYZ(:,1)-       (0))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j) &wholeNodesXYZ(:,2)== (5-0.215/2));
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,4) = ...
              find(abs(wholeNodesXYZ(:,1)-(-0.215/2))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j) &wholeNodesXYZ(:,2)== (5-0.215/2));
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,5) = ...
              i-1 + (j)*3*(InterWallNoElem-1) + (InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,6) = ...
              i-1 + (j)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,7) = ...
              find(abs(wholeNodesXYZ(:,1)-       (0))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j+1) &wholeNodesXYZ(:,2)== (5-0.215/2));
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,8) = ...
              find(abs(wholeNodesXYZ(:,1)-(-0.215/2))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j+1) &wholeNodesXYZ(:,2)== (5-0.215/2));
          
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,1) = ...
              i-1 + (j-1)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,2) = ...
              i-1 + (j-1)*3*(InterWallNoElem-1) + (InterWallNoElem-1)*2 + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,3) = ...
              find(abs(wholeNodesXYZ(:,1)-(+0.215/2))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j) &wholeNodesXYZ(:,2)== (5-0.215/2));
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,4) = ...
              find(abs(wholeNodesXYZ(:,1)-       (0))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j) &wholeNodesXYZ(:,2)== (5-0.215/2));
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,5) = ...
              i-1 + (j)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,6) = ...
              i-1 + (j)*3*(InterWallNoElem-1) + (InterWallNoElem-1)*2 + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,7) = ...
              find(abs(wholeNodesXYZ(:,1)-(+0.215/2))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j+1) &wholeNodesXYZ(:,2)== (5-0.215/2));
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,8) = ...
              find(abs(wholeNodesXYZ(:,1)-       (0))<1e-10 &wholeNodesXYZ(:,3)== heightZ(j+1) &wholeNodesXYZ(:,2)== (5-0.215/2)); 
           
       else
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,1) = ...
              i-1 + (j-1)*3*(InterWallNoElem-1) + (InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,2) = ...
              i-1 + (j-1)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,3) = ...
              i   + (j-1)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,4) = ...
              i   + (j-1)*3*(InterWallNoElem-1) + (InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,5) = ...
              i-1 + (j)*3*(InterWallNoElem-1) + (InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,6) = ...
              i-1 + (j)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,7) = ...
              i   + (j)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(1+(i-1)*2 + (j-1)*InterWallNoElem*2,8) = ...
              i   + (j)*3*(InterWallNoElem-1) + (InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,1) = ...
              i-1 + (j-1)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,2) = ...
              i-1 + (j-1)*3*(InterWallNoElem-1) + (InterWallNoElem-1)*2 + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,3) = ...
              i   + (j-1)*3*(InterWallNoElem-1) + (InterWallNoElem-1)*2 + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,4) = ...
              i   + (j-1)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,5) = ...
              i-1 + (j)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,6) = ...
              i-1 + (j)*3*(InterWallNoElem-1) + (InterWallNoElem-1)*2 + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,7) = ...
              i   + (j)*3*(InterWallNoElem-1) + (InterWallNoElem-1)*2 + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);
          interElem2n(2+(i-1)*2 + (j-1)*InterWallNoElem*2,8) = ...
              i   + (j)*3*(InterWallNoElem-1) + size(wholeNodesXYZ,1)-size(InterNodesXYZ,1);           
       end
    end
end

%% Move the mesh down by 0.5 meter
wholeNodesXYZ(:,3) = wholeNodesXYZ(:,3) - 0.5;
PlotMesh(wholeNodesXYZ, interElem2n, 0)
wholeElem2n = [wholeElem2n; interElem2n];
PlotMesh(wholeNodesXYZ, wholeElem2n, 0)
bottom_elem = [];
for i = 1:length(wholeElem2n)
    if wholeNodesXYZ(wholeElem2n(i,1),3) == min(wholeNodesXYZ(:,3))
        bottom_elem = [bottom_elem,i];
    end
end
PlotMesh(wholeNodesXYZ(:,1:2), wholeElem2n(bottom_elem,1:4),0)

%% InterFace elems
% elem2nInter = elem2nFoundReorder(end-(WallThickNoElem+foundationOutNoElem+foundationInNoElem)*...
%     sum(WallNoElem)+1:end,:);
% PlotMesh(wholeNodesXYZ, elem2nInter, 0)

bottom_nodes = find(wholeNodesXYZ(:,3) == min(wholeNodesXYZ(:,3)));
bottomNodesXYZ = wholeNodesXYZ(bottom_nodes,:);
outerDummyNodesXYZ = offsetWallNodes(wallNodesXYZ(1:sum(WallNoElem),:),...
    wallNodesXYZ(1+sum(WallNoElem):sum(WallNoElem)*2,:),...
    WallThick(1)/WallThickNoElem,...
    (foundationWidth(1)-WallThick(1))/2 +(foundationWidth(1)-WallThick(1))/2/foundationOutNoElem);
innerDummyNodesXYZ = offsetWallNodes(wallNodesXYZ(1+sum(WallNoElem)*WallThickNoElem:sum(WallNoElem)+sum(WallNoElem)*WallThickNoElem,:),...
        wallNodesXYZ(1+sum(WallNoElem)*(WallThickNoElem-1):sum(WallNoElem)*WallThickNoElem,:),...
        WallThick(1)/WallThickNoElem,...
        (foundationWidth(1)-WallThick(1))/2+(foundationWidth(1)-WallThick(1))/2/foundationInNoElem);
interNodesXYZ = [outerDummyNodesXYZ;bottomNodesXYZ;innerDummyNodesXYZ];

elem2nInter = zeros((WallThickNoElem+foundationOutNoElem+foundationInNoElem+2)*...
    sum(WallNoElem),4);
for j= 1:foundationOutNoElem+WallThickNoElem+foundationInNoElem+2
    for i = 1:sum(WallNoElem)
       elem2nInter(i+(j-1)*sum(WallNoElem),1) = i+(j-1)*sum(WallNoElem);
       elem2nInter(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-1)*sum(WallNoElem);
       elem2nInter(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       elem2nInter(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       if i==sum(WallNoElem)
           elem2nInter(i+(j-1)*sum(WallNoElem),1) = i + (j-1)*sum(WallNoElem);
           elem2nInter(i+(j-1)*sum(WallNoElem),2) = 1 + (j-1)*sum(WallNoElem);
           elem2nInter(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem);
           elem2nInter(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem);
       end  
    end
end

% outerDummyCornersXYZ = offsetCurveOut(baseCornersXYZ, baseWall, (foundationWidth-WallThick)/2+(foundationWidth-WallThick)/2/foundationOutNoElem);
% outerDummyNodesXYZ = meshOnCurve(outerDummyCornersXYZ, baseWall, WallNoElem);
% innerDummyCornersXYZ = offsetCurveIn(baseCornersXYZ, baseWall, (foundationWidth-WallThick)/2 + WallThick +(foundationWidth-WallThick)/2/foundationInNoElem);
% innerDummyNodesXYZ = meshOnCurve(innerDummyCornersXYZ, baseWall, WallNoElem);
% interNodesXYZ = [outerDummyNodesXYZ;bottomNodesXYZ;innerDummyNodesXYZ];

 
% elem2nInter = zeros((WallThickNoElem+foundationOutNoElem+foundationInNoElem+2)*...
%     sum(WallNoElem),4);
% for j= 1:foundationOutNoElem+WallThickNoElem+foundationInNoElem+2
%     for i = 1:sum(WallNoElem)
%        elem2nInter(i+(j-1)*sum(WallNoElem),1) = i+(j-1)*sum(WallNoElem);
%        elem2nInter(i+(j-1)*sum(WallNoElem),2) = i+1 + (j-1)*sum(WallNoElem);
%        elem2nInter(i+(j-1)*sum(WallNoElem),3) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem);
%        elem2nInter(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem);
% %        elem2nFound(i+(j-1)*sum(WallNoElem),5) = i+(j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
% %            + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
% %        elem2nFound(i+(j-1)*sum(WallNoElem),6) = i+1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
% %            + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
% %        elem2nFound(i+(j-1)*sum(WallNoElem),7) = i+1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
% %            + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
% %        elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
% %            + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
% %        if j==foundationOutNoElem
% %            elem2nFound(i+(j-1)*sum(WallNoElem),3) = i+1;
% %            elem2nFound(i+(j-1)*sum(WallNoElem),4) = i;
% %            elem2nFound(i+(j-1)*sum(WallNoElem),7) = i+1 + size(wallNodesXYZ,1)...
% %                + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
% %                + sum(WallNoElem)*foundationOutNoElem;
% %            elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + size(wallNodesXYZ,1)...
% %                + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
% %                + sum(WallNoElem)*foundationOutNoElem;
% %        end
%        if i==sum(WallNoElem)
%            elem2nInter(i+(j-1)*sum(WallNoElem),1) = i + (j-1)*sum(WallNoElem);
%            elem2nInter(i+(j-1)*sum(WallNoElem),2) = 1 + (j-1)*sum(WallNoElem);
%            elem2nInter(i+(j-1)*sum(WallNoElem),3) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem);
%            elem2nInter(i+(j-1)*sum(WallNoElem),4) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem);
% %            elem2nFound(i+(j-1)*sum(WallNoElem),5) = i + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
% %                + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
% %            elem2nFound(i+(j-1)*sum(WallNoElem),6) = 1 + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
% %                + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
% %            elem2nFound(i+(j-1)*sum(WallNoElem),7) = 1 + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
% %                + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
% %            elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + sum(WallNoElem) + (j-1)*sum(WallNoElem)+ size(wallNodesXYZ,1)...
% %                + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem);
% %            if j==foundationOutNoElem
% %                elem2nFound(i+(j-1)*sum(WallNoElem),3) = 1;
% %                elem2nFound(i+(j-1)*sum(WallNoElem),4) = i;
% %                elem2nFound(i+(j-1)*sum(WallNoElem),7) = 1 + size(wallNodesXYZ,1)...
% %                    + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
% %                    + sum(WallNoElem)*foundationOutNoElem;
% %                elem2nFound(i+(j-1)*sum(WallNoElem),8) = i + size(wallNodesXYZ,1)...
% %                    + sum(WallNoElem)*(foundationInNoElem+foundationOutNoElem)...
% %                    + sum(WallNoElem)*foundationOutNoElem;
% %            end
%        end  
%     end
% end
PlotMesh(interNodesXYZ(:,1:2), elem2nInter, 0)
interNodesXYZ(:,3) = 0;
%

% %% --- OPTIONAL END-ROTATION ---
% if applyRotation
%     wholeNodesXYZ = wholeNodesXYZ * Rz.';      % all 3D nodes
%     if exist('interNodesXYZ','var') && ~isempty(interNodesXYZ)
%         interNodesXYZ = interNodesXYZ * Rz.';  % interface nodes
%     end
%     baseCornersXYZ = baseCornersXYZ * Rz.';    % reference corners for post-processing
% end

%% quick visual check

PlotMesh(wholeNodesXYZ, wholeElem2n, 0);
axis equal; grid on; view(3);
title(sprintf('AFTER rotation: %.1f° about Z', theta));

%% Save
%%oldFolder = cd(oldFolder);
% timeString = datestr(datetime('now'),30);
% timberElemIndex = [1 1];
% title = strcat('burd_mat_rotated_',timeString,'.mat');
% save(title,'wholeNodesXYZ','wholeElem2n','interNodesXYZ','elem2nInter', 'timberElemIndex')