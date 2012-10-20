#include <iostream>
#include <cstdlib>
#include <cmath>
#include <boost/math/constants/constants.hpp>

using namespace std ;

main()
{
  double coefft ;
  int n,c;
  double pi = boost::math::constants::pi<double>();
  for (n=1; n<=100; n++) {
    /* sawtooth */
    /*if (n%2 == 1) {
      c = n / 2 + 1;
      coefft = -1 / (c * pi) ;
    } else {
      coefft = 0 ;
    }
    cout << coefft << endl ;*/
    /* square */
    if (n%4 == 1) {
      c = n /2 + 1;
      coefft = 4 / (c * pi);
    } else {
      coefft = 0 ;
    }
    cout << coefft << endl ;
  }
}

