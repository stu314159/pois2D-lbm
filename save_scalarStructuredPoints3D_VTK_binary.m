function save_scalarStructuredPoints3D_VTK_binary(filename,dataname,...
    data_set,origin,spacing)

[nx,ny,nz]=size(data_set);

% open the file
fid = fopen(filename,'w');

% ASCII file header
fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'VTK from Matlab\n');
fprintf(fid,'BINARY\n\n');
fprintf(fid,'DATASET STRUCTURED_POINTS\n');
fprintf(fid,'DIMENSIONS %d %d %d \n',nx,ny,nz);
fprintf(fid,'ORIGIN  %4.3f   %4.3f  %4.3f \n',origin(1),origin(2),origin(3));
fprintf(fid,'SPACING %4.3f   %4.3f  %4.3f \n',spacing(1),spacing(2),spacing(3));
fprintf(fid,'\n');
fprintf(fid,'POINT_DATA %d \n',nx*ny*nz);
fprintf(fid,strcat('SCALARS','\t ',dataname,' float ', '\n'));
%fprintf(fid,'SCALARS displacement float \n');
fprintf(fid,'LOOKUP_TABLE default \n');
%fprintf(fid,'\n');

fwrite(fid,reshape(data_set,1,nx*ny*nz),'float','b');
fclose(fid);

return
