#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <cmath>

#define eps std::numeric_limits<double>::epsilon()

using namespace std;

double sawtooth(double x,double timePeriod,double peak,double slope)
{
  double fx ;
  if (x > 0 && x < timePeriod)
  {
    fx = x * peak / timePeriod ;
    //cout << x << "\t" << fx << endl;
    return fx;
  }
  else if (x > timePeriod)
  {
    while (x > timePeriod)
      x = x - timePeriod ;
    return sawtooth(x,timePeriod,peak,slope) ;
  }
  else if (x < 0)
  {
    while (x < 0)
    {
      if (fabs(x) > eps)
        x = x + timePeriod ;
      else 
        return x * peak / timePeriod ;
    }
    return sawtooth(x,timePeriod,peak,slope) ;
  }
  else return 0 ;
}

main()
{
  cout << "eps = " << eps << endl ;
  double T = 0.1;
  double peak = 1.0;
  double slope = peak/T; 

  double x,y ;
  ofstream saw;
  saw.open("fun.txt");
  for (x=-1;x<=1;x+=0.03)
  {
    y = sawtooth(x,T,peak,slope) ;
    saw << x << "\t" << y << endl;
  }
  saw.close();
  system("gnuplot -persist plotScript.p");
}
