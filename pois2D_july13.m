%pois2D_july12.m

% for Pouiseuille flow, I want to be able to measure

% 1. Numerical convergence
% 2. Stabilization rate

% compare both for Zou/He and Regularized BCs
% compare both for LBGK and MRT 


% provide parabolic inlet velocity profile for convergence testing.
% provide uniform inlet velocity profile to show flow development.


clear
clc
close('all')

location = 'office';
% 'home', 'office', 'hamming'

dynamics = 1;
% 1 =LBGK
% 2 = RBGK
% 3 = MRT

Num_ts = 30000;
ts_rep_freq = 1000;
plot_freq = 5000;

Re = 10;
dt = 1e-3;
Ny_divs = 30;

obst_type = 'none';
% 'none'

sim_name = 'pois2D_convergence';
ts_num=0;

Lx_p = 10;
Ly_p = 1;

switch obst_type 
    
    case 'none'
        
        Lo = Ly_p;
        
end


fluid = 4;

switch fluid
    case 1
        rho_p = 1260;
        nu_p = 1.49/rho_p;
        
    case 2
        rho_p = 965.3;
        nu_p = 0.06/rho_p;
        
    case 3
        rho_p = 1000;
        nu_p = 1e-3/rho_p;
        
    case 4
        rho_p = 1000;
        nu_p = 0.001;
        
end

% non-dimensionalize and set up LBM lattice
Uo = nu_p*Re/Lo;
To = Lo/Uo;
Uavg = Uo;

Ld = 1; Td = 1; Ud = (To/Lo)*Uavg;
nu_d = 1/Re;

dx = 1/(Ny_divs-1);
u_lbm = (dt/dx)*Ud;
nu_lbm=(dt/(dx^2))*nu_d;
omega = get_BGK_Omega(nu_lbm);

u_conv_fact = (dt/dx)*(To/Lo);
t_conv_fact = (dt*To);
l_conv_fact = dx*Lo;
p_conv_fact = ((l_conv_fact/t_conv_fact)^2)*(1/3); % <--for EOS type methods...

rho_lbm = rho_p;
rho_out = rho_lbm;

% generate LBM lattice
xm = 0; xp = Lx_p;
ym = 0; yp = Ly_p;

Ny = ceil((Ny_divs-1)*(Ly_p/Lo))+1;
Nx = ceil((Ny_divs-1)*(Lx_p/Lo))+1;

[gcoord,~,~]=RecMesh(xm,xp,ym,yp,Nx,Ny);
[nnodes,~]=size(gcoord);

x_space = linspace(xm,xp,Nx);
y_space = linspace(ym,yp,Ny);

[X,Y]=meshgrid(x_space,y_space);

[w,ex,ey,bb_spd]=D2Q9_lattice_parameters();
%stream_tgt = genTargetVecD2Q9r2(Nx,Ny);
LatticeSize = [Nx Ny];
LatticeSpeeds = [ex; ey];
stm = genStreamTgtMat(LatticeSize,LatticeSpeeds);

numSpd=9;
M = getMomentMatrix('D2Q9');

switch dynamics
    
    case 1
        S = omega.* eye(numSpd);
    case 2
        % TRT model as described in Kevin Tubbs' dissertation section 3.4
        S = zeros(numSpd);
        S(2,2)= omega;
        S(3,3)=omega;
        S(8,8)=omega;
        S(9,9)=omega;
        
        t_s = (1/2)+1/(12*((1/omega)-0.5));
       
        S(5,5)=1/t_s;
        S(7,7)=1/t_s;
    case 3
        % parameters taken from 
        % Chinese Physics Vol 15 No 8 Aug 2006
        % Simulating high Reynolds number flow in 2D lid driven cavity by
        % MRT etc...
        S = zeros(numSpd);
        S(2,2)=1.1;
        S(3,3)=1.0;
        S(5,5)=1.2;
        S(7,7)=1.2;
        S(8,8)=omega;
        S(9,9)=omega;
        
        
end

omega_op = M\(S*M);

% set node lists 

snl=find((gcoord(:,2)==ym) | (gcoord(:,2)==yp));
inl=find(gcoord(:,1)==xm);
inl=setxor(inl,intersect(snl,inl)); % eliminate solid nodes from inlet list
onl=find(gcoord(:,1)==xp);
onl=setxor(onl,intersect(snl,onl)); % eliminate solid nodes from outlet list

% modify for obstructions
switch obst_type
   
    case 'none'
        
        % nothing to do here...
        
 
    
end


Umax = (3/2)*u_lbm;
by = Ly_p/2;
ux_p_in = Umax*(1-((gcoord(inl,2)-by)/by).^2);
uy_p_in = zeros(length(ux_p_in),1);

ux_theory = [0;ux_p_in;0]; ux_theory = (ux_theory')./u_conv_fact;

pltX = ceil(3*Nx/4);
pltX_xcoord= x_space(pltX);
pltX_list = find(gcoord(:,1)==pltX_xcoord);

[fIn,fOut,rho,ux,uy]=Initialize_F_zero(gcoord,ex,ey,w,rho_lbm);

fEq = zeros(nnodes,numSpd);

v_data = zeros(Num_ts,1);
v_data_LP = ceil(Ny/2)*Nx+ceil(Nx/2);

u_data = zeros(Num_ts,1);
p_ref = rho_out*p_conv_fact;


% prepare for regularized BCs and dynamics
e_i = [ex;ey];

Q_mn = zeros(2,2,9);
for i = 1:9
   Q_mn(:,:,i)=e_i(:,i)*e_i(:,i)' - (1/3)*eye(2,2); 
end

Q_flat = zeros(9,4);
for i = 1:9
   q_tmp = Q_mn(:,:,i); q_tmp = q_tmp(:); q_tmp = q_tmp';
   Q_flat(i,:) = q_tmp;
end


indir_p = find(e_i(1,:)==-1);
indir_0 = find(e_i(1,:)==0);
indir_m = find(e_i(1,:)==1);

outdir_p = find(e_i(1,:)==1);
outdir_0 = find(e_i(1,:)==0);
outdir_m = find(e_i(1,:)==-1);



fprintf('Number of Lattice-points = %d.\n',nnodes);
fprintf('Number of time-steps = %d. \n',Num_ts);
fprintf('LBM viscosity = %g. \n',nu_lbm);
fprintf('LBM relaxation parameter (omega) = %g. \n',omega);
fprintf('LBM flow Mach number = %g. \n',u_lbm);

input_string = sprintf('Do you wish to continue? [Y/n] \n');

run_dec = input(input_string,'s');

if ((run_dec ~= 'n') && (run_dec ~= 'N'))
    
    fprintf('Ok! Cross your fingers!! \n');
    
    
    % add paths for Jacket libraries
     switch location
       
        case 'home'
            addpath('/usr/local/jacket/engine');
            addpath('/home/stu/Dropbox/matlab/jacketSDK/pc_pois2D_velBCs');
            addpath('/home/stu/Dropbox/matlab/jacketSDK/bounce_back_jkt');
            addpath('/home/stu/Dropbox/matlab/jacketSDK/stream_jkt');
            addpath('/home/stu/Dropbox/matlab/jacketSDK/pois2D_LBGK_ts');
            addpath('/home/stu/Dropbox/matlab/jacketSDK/channel2D_VW_PE_LBGK_ts');
            
        case 'office'
            addpath('/usr/local/jacket/engine');
            addpath('/home/srblair/Dropbox/matlab/jacketSDK/pc_pois2D_velBCs');
            addpath('/home/srblair/Dropbox/matlab/jacketSDK/bounce_back_jkt');
            addpath('/home/srblair/Dropbox/matlab/jacketSDK/stream_jkt');
            addpath('/home/srblair/Dropbox/matlab/jacketSDK/pois2D_LBGK_ts');
            addpath('/home/srblair/Dropbox/matlab/jacketSDK/channel2D_VW_PE_LBGK_ts');
        case 'hamming'
        
     end
     
     % send data to the GPU
     
%      fIn = gsingle(fIn);
%      fOut = gsingle(fOut);
%      fEq = gsingle(fEq);
%      ux_p_h = zeros(nnodes,1); ux_p_h(inl)=ux_p_in; ux_p_h(onl)=ux_p_in;
%      ux_p = gsingle(ux_p_h);
%      rho = gzeros(nnodes,1);
%      ux = gzeros(nnodes,1);
%      uy = gzeros(nnodes,1);
%      inl_i = zeros(nnodes,1); inl_i(inl)=1; inl_d = gint32(inl_i);
%      onl_i = zeros(nnodes,1); onl_i(onl)=1; onl_d = gint32(onl_i);
%      snl_i = zeros(nnodes,1); snl_i(snl)=1; snl_d = gint32(snl_i);
     
     
     % prep stuff for writing to vtk
     uz_h = zeros(nnodes,1);
     gcoord_z = zeros(nnodes,1);
     u_data = zeros(Num_ts,1);
     v_data_LP = ceil(Ny/2)*Nx+ceil(Nx/2);
     %stm_d = gint32(stm);
     %bb_spd_d = gint32(bb_spd);
     
     %commence time stepping
     tic;
     for ts = 1:Num_ts
         % say something comforting about my progress...
         if(mod(ts,ts_rep_freq)==0)
             fprintf('Executing time step number %d.\n',ts);
         end
         
         % compute density
        rho = sum(fIn,2);
        
        % compute velocities
        ux = (fIn*ex')./rho;
        uy = (fIn*ey')./rho;
        
        % 
        ux(snl)=0; uy(snl)=0;
        
        uy(inl)=0; ux(inl)=ux_p_in;
        rho(inl)=(1./(1-ux(inl))).*(2*sum(fIn(inl,indir_p),2)+sum(fIn(inl,indir_0),2));
        rho(onl)=rho_out;
        uy(onl)=0;
        ux(onl)=-1+((2*sum(fIn(onl,outdir_p),2)+sum(fIn(onl,outdir_0),2))./rho(onl));
        
        % compute fEq
        for i = 1:numSpd
            cu = 3*(ex(i)*ux+ey(i)*uy);
            fEq(:,i)=w(i)*rho.*(1+cu+(1/2)*(cu.*cu) - ...
                (3/2)*(ux.^2 + uy.^2 ));
        end
        
        fIn(inl,indir_m)=fEq(inl,indir_m)+fIn(inl,bb_spd(indir_m))-fEq(inl,bb_spd(indir_m));
        fIn(onl,outdir_m)=fEq(onl,outdir_m)+fIn(onl,bb_spd(outdir_m))-fEq(onl,bb_spd(outdir_m));
        
        % now, compute Pi_1 based on these trial values and correct
        Q_tmp = (fIn(inl,:)-fEq(inl,:))*Q_flat;
        f1 = (Q_tmp*(Q_flat').*repmat(w,length(inl),1))*(9/2);
        fIn(inl,:)=fEq(inl,:)+f1;
        
        Q_tmp = (fIn(onl,:)-fEq(onl,:))*Q_flat;
        f1 = (Q_tmp*(Q_flat').*repmat(w,length(onl),1))*(9/2);
        fIn(onl,:)=fEq(onl,:)+f1;
        
        switch dynamics
            
            case 1
                fOut = fIn - omega*(fIn-fEq);
                
            case 2
                % collision
                Q_tmp=(fIn-fEq)*Q_flat;
                f1 = (Q_tmp*(Q_flat').*repmat(w,nnodes,1))*(9/2);
                f1 = f1+fEq;
                fOut = f1-(f1-fEq)*omega;
                
            case 3
                fOut = fIn - (fIn - fEq)*omega_op;
        end
         
        % bounce-back
        for i = 1:numSpd
            fOut(snl,i)=fIn(snl,bb_spd(i));
        end
        
        
        % stream
        %fIn(stream_tgt)=fOut(:);
        for i = 1:numSpd
            fIn(stm(:,i),i)=fOut(:,i);
        end
        
        
        u_data(ts)=ux(v_data_LP)./u_conv_fact;
        
         if(mod(ts,plot_freq)==0)
             % plot something....
             % plot something, plot something cool!!
            ux_h = (ux./u_conv_fact);
            uy_h = (uy./u_conv_fact);
            %uz_h = (uz./u_conv_fact);
            
            pressure_h = ((rho)*p_conv_fact-p_ref);
            
            vtk_suffix=sprintf('_velocityAndPressure%d.vtk',ts_num);
            ts_fileName=strcat(sim_name,vtk_suffix);
            save_velocityAndPressureVTK_binaryR2(pressure_h,ux_h,uy_h,uz_h,...
                 gcoord(:,1),gcoord(:,2),gcoord_z,ts_fileName,[Nx Ny 1]);
            
            ts_num=ts_num+1;
             
         end
         
     end
     
     ex_time = toc;
     fprintf('LPU/sec = %g.\n',Num_ts*nnodes/ex_time);
     
     % plot velocity profile and compare against theory
     ux_lbm = ux_h(pltX_list);
     figure(1)
     plot(ux_theory,y_space,'-r',ux_h(pltX_list),gcoord(pltX_list,2),'xb');
     
     rel_err = norm(ux_theory - ux_lbm',2)/norm(ux_theory,2);
     fprintf('Relative error = %g.\n',rel_err);
     fprintf('Grid Resolution = %d.\n',Ny_divs);
     fprintf('dx = %g.\n',dx);
     fprintf('dt = %g.\n',dt);
     fprintf('omega = %g.\n',omega);
     
     figure(2)
     plot(1:Num_ts,u_data,'LineWidth',2);
     grid on
     title('\bf{Horizontal Velocity Fluctuation x/Lx = 0.75 vs Time Step, Re = 10}','FontSize',12);
     xlabel('\bf{Time Step}','FontSize',12);
     ylabel('\bf{Umax (m/sec)}','FontSize',12);
else
    fprintf('Run aborted.  Better luck next time!\n');
end  