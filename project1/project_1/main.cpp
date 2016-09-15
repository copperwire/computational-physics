#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <iomanip>

using namespace std;

double funct(double x){
    return 100*exp(-10 * x);
}
double exact(double x){
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}


void non_special_method(double U [], int exponent){
    char filename []  = "forw_back.txt";
    int N = pow(10, exponent);
    double h = 1./(N+1);
    double hh = h*h;

    //A = upper diagonal
    //B = lower diagonal
    //M = midle diagonal
    //F = function evaluation
    //U = unknown function result

    double* A = new double[N+2];
    double* B = new double[N+2];
    double* M = new double[N+2];
    double* F = new double[N+2];



    /* Remember to write to file
     * Implement the exact function in python - faster
     * NB!
     * The function needs to be EXACTLY 0 at boundary conditions DIRICHLET MOTHERFUCKER
    c*/

    ofstream datafile (filename);

    for(int i = 0; i<N+2; i++){
        A[i] = -1.0;
        B[i] = -1.0;
        M[i] = 2.0;
        F[i] = hh*100*exp(-10* i*h);

    }

    // Equations are fetched from the lecture notes on linear algebra regard. a tri-diagonal matrix on a general form.

    for(int i = 2; i < N +2; i++){
        M[i] = M[i] - (B[i] * A[i-1])/M[i-1];
        F[i] = F[i] - (B[i] * F[i-1])/M[i-1];
    }

    // writing outside loop to include boundary conditions
    datafile << U[0] << endl;

    for(int i = N + 1; i > 1; i--){
        U[i - 1] = (F[i-1] -  A[i-1] * U[i])/M[i-1] ;
        datafile << U[i -1] << endl ;
    }

    datafile << U[N+2] ;
}


void special_method(double U [], int exponent){
    for(int i = 1; i <= exponent; i++){
        string filename  = "special_result";
        string rel_err_file = "rel_err";

        string new_string = filename + rel_err_file;

        string argument = to_string(i);
        filename.append(argument);
        rel_err_file.append(argument);

        int N = pow(10, i);
        double h = 1./(N+1);
        double hh = h*h;

        //M = midle diagonal
        //F = function evaluation
        //U = unknown function result

        double* M = new double[N+2];
        double* F = new double[N+2];
        double* x = new double[N+2];


        M[1] = M[N] =  2.0;
        F[1] = hh*100*exp(-10 * h);

        for (int i = 2; i< N+2; i++){
            x[i] = i*h;
        }


        for(int i = 2; i < N+2 ; i++){
            M[i] = (i + 1.0)/((double) i);
            F[i] = hh*funct(x[i]) + F[i-1]/M[i-1] ;
        }

         ofstream datafile (filename.c_str());

         datafile << U[N+1] << endl;

         U[N] = F[N]/M[N];

         datafile << U[N] << endl;

         for(int i = N-1 ; i > 0; i--){
             U[i] = (F[i] + U[i+1])/M[i];
             datafile << U[i] << endl ;
         }

         datafile << U[0] ;


         ofstream rel_file (rel_err_file.c_str());

         for(int i = 2; i < N+2; i++){
            double exact_val  = exact(x[i]);
            double Relerr = fabs( (exact_val - U[i]) / exact_val);

            rel_file << setw(15) << setprecision(8) << x[i];
            rel_file << setw(15) << setprecision(8) << U[i];
            rel_file << setw(15) << setprecision(8) << exact_val;
            rel_file << setw(15) << setprecision(8) << log10(Relerr) << endl;
         }


    }

}


int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    double* U = new double[N+2];
    double* U_spec = new double[N+2];

    non_special_method(U, N);
    special_method(U_spec, N);

    return 0;

}
