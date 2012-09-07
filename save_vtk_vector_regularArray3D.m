function save_vtk_vector_regularArray3D(X1,X2,X3,filename,vector_name,origin,spacing)

[nx,ny,nz]= size(X1);
fid = fopen(filename,'wt'); %open the file in text mode
fprintf(fid,'# vtk DataFile Version 2.0 \n');
fprintf(fid,'Comment goes here...\n');
fprintf(fid,'ASCII \n');
fprintf(fid,'\n');
fprintf(fid,'DATASET STRUCTURED_POINTS\n');
fprintf(fid,'DIMENSIONS %d   %d  %d \n',nx,ny,nz);
fprintf(fid,'\n');
fprintf(fid,'ORIGIN   %4.3f    %4.3f   %4.3f \n',origin(1),origin(2),origin(3));
fprintf(fid,'SPACING  %4.3f    %4.3f   %4.3f \n',spacing(1),spacing(2),spacing(3));
fprintf(fid, '\n');
fprintf(fid,'POINT_DATA   %d \n', nx*ny*nz);

fprintf(fid,strcat('VECTORS','\t  ',vector_name,' double ','\n'));

fprintf(fid, '\n');
for a=1:nz
    for b=1:ny
        for c=1:nx
            fprintf(fid, '%f ', X1(c,b,a));
            fprintf(fid, '%f ', X2(c,b,a));
            fprintf(fid, '%f ', X3(c,b,a));
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);
return