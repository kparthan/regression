#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std ;

main()
{
	ifstream file ("samples.txt") ;
	int a[10] ;
	for (int i=0; i<10; i++)
		file >> a[i] ;
	for (int i=0; i<10; i++)
		cout << a[i] << endl ;
}

