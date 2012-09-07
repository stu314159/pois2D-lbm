function save_vectorStructuredGrid3D_VTK_binary(filename,...
    gcoord,d_x,d_y,d_z,data_shape)

[numEl,~]=size(gcoord);
%[nx,ny,nz]=size(data_set);


fid = fopen(filename,'w');

%ASCII file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'BINARY\n\n');
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid,['DIMENSIONS ' num2str(data_shape(1)) ' ' num2str(data_shape(2)) ' ' num2str(data_shape(3)) '\n']);
fprintf(fid,['POINTS ' num2str(numEl) ' float\n']);
fclose(fid);
% append binary xyz data
fid = fopen(filename,'a');
fwrite(fid,[gcoord(:,1)'; gcoord(:,2)'; gcoord(:,3)'],'float','b');

% append vector data
fwrite(fid,[reshape(d_x,1,numEl); reshape(d_y,1,numEl); reshape(d_z,1,numEl)],'float','b');

fclose(fid);