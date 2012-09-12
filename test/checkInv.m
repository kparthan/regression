phi = load('../temp/phi') ;
size(phi) ;

phiT = phi' ;
phiTphi = phiT * phi ;
det(phiTphi) 

inverse = inv(phiTphi) ;
I1 = phiTphi * inverse ;
printMatrix(I1,'../temp/I1') ;

invFile = load('../temp/inverse') ;
I2 = phiTphi * invFile ;
printMatrix(I2,'../temp/I2') ;

diff = inverse - invFile ;
printMatrix(diff,'../temp/invDiff') ;

