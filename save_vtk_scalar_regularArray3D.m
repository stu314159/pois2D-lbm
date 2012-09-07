function save_vtk_scalar_regularArray3D(array,filename,scalar_name,origin,spacing)

[nx,ny,nz]= size(array);
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
fprintf(fid,strcat('SCALARS','\t  ',scalar_name,' double ','\n'));
%fprintf(fid,'SCALARS pressure double \n');
fprintf(fid,'LOOKUP_TABLE default \n');
fprintf(fid, '\n');
for a=1:nz
    for b=1:ny
        for c=1:nx
            fprintf(fid, '%d ', array(c,b,a));
        end
        fprintf(fid, '\n');
    end
end


fclose(fid);  % close the file like a good boy

return