%phi = load('phi_ublas') ;
phi = load('phi') ;
disp('size(phi): ') ;
size(phi) 

phi_trans = phi' ;
A = phi_trans * phi ;
A_inv = inv(A) ;
disp('size(A_inv): ') ;
size(A_inv) 

%inv_computed = load('inverse_ublas') ;
inv_computed = load('inverse') ;
inv_diff = inv_computed - A_inv ;

%fd = fopen('check/diff_inv_ublas_50.txt','w') ;
fd = fopen('check/diff_inv_50.txt','w') ;

for i=1:size(A_inv,1)
    for j=1:size(A_inv,2)
        fprintf(fd,'%f ',inv_diff(i,j)) ;
    end    
    fprintf(fd,'\n') ;
end    

fclose(fd) ;
