function [] = print (fileName,matrix)

    fw = fopen(fileName,'w');
    for i=1:size(matrix,1)
       for j=1:size(matrix,2)
          fprintf(fw,'%f ',matrix(i,j)) ; 
       end    
       fprintf(fw,'\n') ;
    end
    
    fclose(fw);
end

