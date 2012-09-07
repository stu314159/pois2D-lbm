function [fIn,fOut,rho,ux,uy] = Initialize_F_Pois(gcoord,L,u_bc,ex,ey,w,rho_init)

[~,numSpd]=size(ex);
[nnodes,~]=size(gcoord);
b = L/2;

fIn = zeros(nnodes,numSpd);
%fOut = fIn;
rho = rho_init.*ones(nnodes,1);

ux = (u_bc)*((1-((gcoord(:,2)-b).^2)./(b*b)));
uy = zeros(nnodes,1);

for i = 1:numSpd
   cu = 3*(ex(i)*ux + ey(i)*uy);
   fIn(:,i) = w(i)*rho.*(1+cu+(1/2)*(cu.*cu)-(3/2)*(ux.^2 + uy.^2));
    
end

fOut = fIn; % just to start 