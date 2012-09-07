function fIn_b = velocityBC_D2Q9(fIn_b,w,ex,ey,ux_p,uy_p)

rho = sum(fIn_b,2);
ux = (fIn_b*ex')./rho;
uy = (fIn_b*ey')./rho;


dx = ux_p-ux;
dy = uy_p-uy;



for spd=2:9
   cu=(3)*(ex(spd)*dx+ey(spd)*dy);
   fIn_b(:,spd)=fIn_b(:,spd)+w(spd)*(rho.*cu);
       
end

