function [] = printMatrix(M,file)

    fw = fopen(file,'w') ;

    rows = size(M,1) ;
    cols = size(M,2) ;

    for i=1:rows
        for j=1:cols
            fprintf(fw,'%f ',M(i,j)) ;
        end
        fprintf(fw,'\n') ;
    end    

    fclose(fw) ;

end