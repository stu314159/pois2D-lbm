% written by: Stu Blair
% date: Jan 16, 2012
% purpose: get force vector for all solid nodes in an LBM domain

function F = getPressureForce2D(solidNodes,streamTgtMat,LatticeSpeeds,bb_spd,fIn,...
    time_conv,space_conv)

% determine number of solid nodes
numSolidNodes=length(solidNodes);

%deterimne number of dimensions
[numDim,numSpd]=size(LatticeSpeeds);

F = zeros(numSolidNodes,numDim);
% get forceContribMat - map of force contributors for each solid node
forceContribMat = zeros(numSolidNodes,numSpd);
for spd = 2:numSpd % stationary speed never contributes
    
    % look to see if neighbor is a solid node or not if it is not, then it
    % contributes its momentum to the force...
    forceContribMat(:,spd)=~ismember(streamTgtMat(solidNodes,spd),solidNodes);
end

e_alpha = LatticeSpeeds';

for spd=2:numSpd
    % sum of density distribution heading out from a particle in a given direction
    %, plus the
    % amount headed in from ngb in corresponding direction
    f1 = fIn(solidNodes,spd)+fIn(streamTgtMat(solidNodes,spd),bb_spd(spd));
    % make the result zero if solid node.
    f2 = (f1*((space_conv/time_conv)^2)*(1/3)*space_conv).*forceContribMat(:,spd);
    % distribute in the appropriate directions...
    f3 = kron(f2,e_alpha(bb_spd(spd),:));
    % accumulate the result
    F = F + f3;
  
end
