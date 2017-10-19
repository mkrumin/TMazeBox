function newMaps = remapROI(maps, source, dest)

% for debugging
% maps = {[1000:10000]};
% source = [10, 20, 480, 500];
% dest = [-50, -150, 450, 400];

nMaps = length(maps); % maps is a cell array

topY = min(source(1), dest(1));
leftX = min(source(2), dest(2));
bottomY = max(source(3), dest(3));
rightX = max(source(4), dest(4));


hSource = source(3)-source(1) + 1;
wSource = source(4)-source(2) + 1;

newMaps = cell(size(maps));

for iMap = 1:nMaps
    oldMask = zeros(hSource, wSource);
    oldMask(maps{iMap}) = 1;
    
    % this mask will include both the source and the destination masks
    largeMask = zeros(bottomY-topY+1, rightX-leftX+1);
    
    idxY = source(1) - topY +1 : source(3) - topY +1;
    idxX = source(2) - leftX +1 : source(4) - leftX +1;
    largeMask(idxY, idxX) = oldMask;
    
    idxY = dest(1) - topY +1 : dest(3) - topY +1;
    idxX = dest(2) - leftX +1 : dest(4) - leftX +1;
    newMask = largeMask(idxY, idxX);
    newMaps{iMap} = find(newMask(:));
end