%%Exportation matrice carrée de MatLab à Mathematica
function toto=export_matrix_MatLab2Mathematica(fichier,matrix)

fid = fopen(fichier,'w');
fprintf(fid,'mat={');
for i=1:(size(matrix,1))
    fprintf(fid,'{');
    for j=1:(size(matrix,1))
        %matrix(i,j)
        fprintf(fid,'%f',matrix(i,j));
        if(j<size(matrix,1))
            fprintf(fid,',');
        end
    end
    fprintf(fid,'}');
    if(i<size(matrix,1))
            fprintf(fid,',');
    end
end
fprintf(fid,'}');

status = fclose(fid);
