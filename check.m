phi = load('temp/phi') ;
y = load('temp/y') ;
disp('size(phi): ') ;
size(phi) ;
format long;
phi_trans = phi' ;
A = phi_trans * phi ;
matlab_inv = inv(A) ;

I2 = A * matlab_inv;
print('temp/I2',I2);

print('temp/matlab_inv',matlab_inv);
disp('size(A_inv): ') ;
size(matlab_inv) ;
weights = matlab_inv * phi_trans * y 

computed_inv = load('temp/inverse') ;

I1 = A * computed_inv;
print('temp/I1',I1);
inv_diff = computed_inv - matlab_inv ;
print('temp/inv_diff',inv_diff);
