function S = getEwMatrixMRT(lattice,omega)

switch lattice
   
    case 'D2Q9'
        S = zeros(9);
       S(2,2)=1.1;
       S(3,3)=1.0;
       S(5,5)=1.2;
       S(7,7)=1.2;
       S(8,8)=omega;
       S(9,9)=omega;
        
       
        
    case 'D3Q15'
        % from D'Humieres et al Phil. Trans. R. Soc. Lond. A (2002) v360,
        % 437-451.
        s1 = 1.6;
        s2 = 1.2;
        s4 = 1.6;
        s14 = 1.2;
        s9 = omega; 
        s11 = omega;
        
        S = diag([0 s1 s2 0 s4 0 s4 0 s4 s9 s9 s11 s11 s11 s14]);
        
    case 'D3Q19'
        s1=1.19;
        s2=1.4; 
        s4=1.2;
        s9=omega;
        s13=omega;
        s16=1.98;
        s10=1.4;
        
        S = diag([0 s1 s2 0 s4 0 s4 0 s4 s9 s10 s9 s10 s13 s13 s13 s16 s16 s16]);
        
        
    otherwise
        disp('invalid lattice type for MRT');
    
    
end