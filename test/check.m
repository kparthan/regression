phi = zeros(100,3);
phi(:,1) = 1.0;

x = 0:0.0628:6.2172;

for i=1:100
    arg = 2 * pi * x(i) ;
    phi(i,2) = sin(arg) ;
    phi(i,3) = cos(arg) ;
end    

fw = fopen('phi_matlab','w') ;
for i=1:100
    for j=1:3
        fprintf(fw,'%f \n',phi(i,j)) ;
    end
    fprintf(fw,'\n') ;
end    

fclose(fw) ;

phiT = phi' ;
phiTphi = phiT * phi ;
det(phiTphi) ;

inverse = inv(phiTphi) ;

y = load('../temp/y') ;
w = inverse * phiT * y


