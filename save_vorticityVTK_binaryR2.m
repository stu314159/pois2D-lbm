function save_vorticityVTK_binaryR2(u,v,w,x,y,z,filename,dims)
nr_of_elements=numel(x);
fid = fopen(filename, 'w'); 

%ASCII file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'BINARY\n\n');
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid, ['DIMENSIONS ' num2str(dims(1)) ' ' num2str(dims(2)) ' ' num2str(dims(3)) '\n']);
fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
fclose(fid);

%append binary x,y,z data
fid = fopen(filename, 'a'); 
fwrite(fid, [reshape(x,1,nr_of_elements);  reshape(y,1,nr_of_elements); reshape(z,1,nr_of_elements)],'float','b');

%append another ASCII sub header
fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
fprintf(fid, 'VECTORS vorticity_vectors float\n');

%append binary u,v,w data
fwrite(fid, [reshape(u,1,nr_of_elements);  reshape(v,1,nr_of_elements); reshape(w,1,nr_of_elements)],'float','b');

%append some scalar data
fprintf(fid, '\nSCALARS VorticityMagnitude float\n'); %ASCII header
fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
fwrite (fid, reshape(sqrt(u.^2+v.^2+w.^2),1,nr_of_elements),'float','b'); %binary data

% %append some scalar data
% fprintf(fid, '\nSCALARS Pressure float\n'); %ASCII header
% fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
% fwrite (fid, reshape(pressure,1,nr_of_elements),'float','b'); %binary data

fclose(fid);


return