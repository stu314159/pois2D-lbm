function S = getEwMatrixMRT(lattice,omega)

switch lattice
   
    case 'D2Q9'
        S = zeros(9);
        S(2,2)=omega;
        S(3,3)=omega;
        S(8,8)=omega;
        S(9,8)=omega;
        
        t_s = (1
        
    case 'D3Q15'
        
    case 'D3Q19'
        
    otherwise
        disp('invalid lattice type for MRT');
    
    
end