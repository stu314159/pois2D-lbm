function stv = genStreamTgtMat(LatticeSize,LatticeSpeeds)


nDims=length(LatticeSize); %format for Lattice size (Nx Ny [Nz])
nnodes = prod(LatticeSize);
ind = 1:nnodes;
[~,numSpd]=size(LatticeSpeeds);
stv = zeros(nnodes,numSpd);

switch nDims
    case 2
        Nx = LatticeSize(1); Ny = LatticeSize(2);
        ind = reshape(ind,[Nx Ny]);
        
        ex = LatticeSpeeds(1,:);
        ey = LatticeSpeeds(2,:);
        
        for spd = 1:numSpd
            t=circshift(ind,[-ex(spd) -ey(spd)]);
            stv(:,spd)=t(:);
        end
         
        
        
    case 3
        Nx = LatticeSize(1); Ny = LatticeSize(2);
        Nz = LatticeSize(3);
        ind =reshape(ind,[Ny Nx Nz]);
        ex = LatticeSpeeds(1,:);
        ey = LatticeSpeeds(2,:);
        ez = LatticeSpeeds(3,:);
        
        for spd = 1:numSpd
            t = circshift(ind,[-ey(spd) -ex(spd) -ez(spd)]);
            stv(:,spd)=t(:);
        end
        
        
end






