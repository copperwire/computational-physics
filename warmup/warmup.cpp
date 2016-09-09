#include <cmath> /* sqrt, atan */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <typeinfo>
using namespace std;


double func_1(double x, double h){
	double result; 
	result = (atan(x + h) - atan(x)) / h ;
	return result;
}

double func_2(double x, double h){
	double result;
	result = (atan(x + h) - atan(x - h))/ (2.0*h) ;
	return result;
}

void write_to_file(ostream datafile, double f1, double f2, double h){

}


int main(int argc, char *argv[])
{
	int N;
	double x, exact, inc, h;
	double* step_vec; 
	double* func1_vec; 
	double* func2_vec;
	char* filename;

	filename =argv[2];
	printf("%s\n", filename);

	x = sqrt(2);
	N = atoi(argv[1]);
	exact = 1.0/3;
	func1_vec = new double[N];
	func2_vec = new double[N];
	step_vec = new double[N];
	inc = 0.00001/N;

	//Need to figure out what type of object filename is?!??!?!
	printf ("interval: %i \n", N);
	printf ("type of filename var %s \n", typeid(filename).name());

	ofstream datafile (filename);	
	for(int i = 1; i <= N+1; i += 1){
		h = inc*i;
		step_vec[i-1] = x;
		func1_vec[i-1] = log(abs((func_1(x, h) - exact)/exact));
		func2_vec[i-1] = log(abs((func_2(x, h) - exact)/exact));
		datafile << func1_vec[i-1] << " ";
		datafile << func2_vec[i-1] << " ";
		datafile << h << " ";
		datafile << "\n";
	}

	datafile.close();
	return 0;
}


