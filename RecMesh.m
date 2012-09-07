% function [GCOORD,NODES,BOUNDARIES]= RecMesh(x_left,x_right,y_bottom,y_top,Nx,Ny)
%
% Written by: Stu Blair
% Date: 2/9/2010
% Purpose: provide mesh definition data structures for a structured,
% uniform rectangular mesh to be used with ME 4612 final project

function [GCOORD,NODES,BOUNDARIES] = RecMesh(x_left,x_right,y_bottom,y_top,Nx,Ny)

x_space = linspace(x_left,x_right,Nx);
y_space = linspace(y_bottom,y_top,Ny);

[X,Y] = meshgrid(x_space,y_space);

nnodes = Nx*Ny;
nnel = 4; %constant for rectangular elements
nel = (Nx-1)*(Ny - 1); 

GCOORD = zeros(nnodes,2);
NODES = zeros(nel,nnel);

% === form the GCOORD array ===============================================
index = 0;
for j = 1:Ny
    for i = 1:Nx
        index = index+1;
        GCOORD(index,1) = X(j,i);
        GCOORD(index,2) = Y(j,i);
    end
end

% =========================================================================

% ====== Form the NODES array =============================================

for j = 1:(Ny - 1)
    for i = 1:(Nx - 1)
        ele_index = (j-1)*(Nx-1)+i;
        node1 = (j-1)*Nx + i;
        node2 = node1+1;
        node3 = node2+Nx;
        node4 = node3-1;
        NODES(ele_index,1) = node1;
        NODES(ele_index,2) = node2;
        NODES(ele_index,3) = node3;
        NODES(ele_index,4) = node4;
    end
end
% =========================================================================

BOUNDARIES.bottom = 1:Nx; %left-to-right
BOUNDARIES.right = Nx:Nx:Nx*(Ny); %bottom-to-top
BOUNDARIES.top = Nx*Ny:-1:(Nx*(Ny-1)+1);%right-to-left
BOUNDARIES.left = (Nx*(Ny-1)+1):-Nx:1;%top-to-bottom

% this array of BOUNDARIES will have ALL of the nodes on the given
% boundary.  This is important for carrying out the boundary integrals on
% each side.  Somehow, in the application of these BCs, I will have to deal
% with the fact that some nodes are shared between boundaries.  For this
% simple geometry, the first and last node on each boundary is shared.

% based on the logic of forming this array, the nodes are ordered such that
% the domain is encircled counter-clockwise starting at the bottom
% left-hand corner.