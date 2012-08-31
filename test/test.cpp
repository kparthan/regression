#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

using namespace std ;

int main(int argc, char **argv)
{
	/*ifstream file ("samples.txt") ;
	int a[10] ;
	for (int i=0; i<10; i++)
		file >> a[i] ;
	for (int i=0; i<10; i++)
		cout << a[i] << endl ;*/
 /* string x ;
	if (x.compare("something") == 0)
		cout << "something" << endl ;
	x = "bjhbbnjknklnmlk" ;
	cout << x << endl ;
	x = "hello" ;
	cout << x << endl ;

	if (string(argv[1]).compare(x) == 0)
	{
		cout << "equal" << endl ;
		strcpy(x,"copied-1") ;
		cout << x << endl ;
	}
	else
	{
		cout << "not equal" << endl ;
		strcpy(x,"copied-2") ;
		cout << x << endl ;
	}*/
	/*int x = 2 + 3 +
					5 + 9 ;
	cout << "x = " << x << endl ;
	cout << "log(2) = " << log2(1024.0) << endl ;*/

  int x[] = {1,4,5,8,77} ;
  vector<int> xv (x, x+sizeof(x)/sizeof(int)) ;
  cout << xv.size() << endl ;
  
  long double y[] = {1.77,4,5,8,7,8798.7} ;
  vector<long double> yv (y, y+sizeof(y)/sizeof(long double)) ;
  cout << yv.size() << endl ;
}

