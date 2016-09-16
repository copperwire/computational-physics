#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <time.h>
#include <armadillo>

using namespace std;

double funct(double x){
    return 100*exp(-10 * x);
}
double exact(double x){
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}

int non_special_method(int exponent){

    string time_file = "general_time";
    ofstream time_write (time_file.c_str());
    time_write << "#num iterations, seconds/clocks_per_sec" <<endl;

    for(int i = 1; i <= exponent; i++){

        string filename = "general_results";
        string argument = to_string(i);
        filename.append(argument);
        filename.append(".txt");

        int N = pow(10, i);
        double h = 1./(N+1);
        double hh = h*h;
        double start, stop;

        //A = upper diagonal
        //B = lower diagonal
        //M = midle diagonal
        //F = function evaluation
        //U = unknown function result

        double* A = new double[N+2];
        double* B = new double[N+2];
        double* M = new double[N+2];
        double* F = new double[N+2];
        double* U = new double[N+2];
        double* x = new double[N+2];

        /* Remember to write to file
         * Implement the exact function in python - faster
         * NB!
         * The function needs to be EXACTLY 0 at boundary conditions DIRICHLET MOTHERFUCKER
        */

        x[1] = U[1] = F[1] = 0; x[N+1] = U[N+1] = F[N+1] = 0;

        for (int i = 1; i<N+2; i++){
            x[i] = i*h;
            F[i] = hh*100*exp(-10* x[i]);
        }
        for(int i = 0; i<N+2; i++){
            A[i] = -1.0;
            B[i] = -1.0;
            M[i] = 2.0;
        }

        // Start of algoritm

        start = clock();
        for(int i = 2; i < N +2; i++){
            M[i] = M[i] - (B[i] * A[i-1])/M[i-1];
            F[i] = F[i] - (B[i] * F[i-1])/M[i-1];
        }

        U[N] = F[N]/M[N];

        for(int i = N -1 ; i > 1; i--){
            U[i] = (F[i] -  A[i] * U[i +1])/M[i] ;
        }
        stop = clock ();
        // end of algorithm

        time_write << N <<"  "<<(double (stop - start) / CLOCKS_PER_SEC) << endl ;

        ofstream rel_file (filename.c_str());

        for(int i = 1; i < N + 2; i++){
           double exact_val  = exact(x[i]);
           double Relerr = fabs( (exact_val - U[i]) / exact_val);

           rel_file << setw(15) << setprecision(8) << x[i];
           rel_file << setw(15) << setprecision(8) << U[i];
           rel_file << setw(15) << setprecision(8) << exact_val;
           rel_file << setw(15) << setprecision(8) << log10(Relerr) << endl;
           }



        rel_file.close();
        delete A; delete B; delete M; delete F; delete U; delete x;
        }
    return 0;
    time_write.close();
}


int special_method(int exponent){
    string time_file = "special_time";
    ofstream time_write (time_file.c_str());
    time_write << "#num iterations, seconds/clocks_per_sec" <<endl;

    for(int i = 1; i <= exponent; i++){

        string rel_err_file = "special_result";

        string argument = to_string(i);
        rel_err_file.append(argument);
        rel_err_file.append(".txt");

        int N = pow(10, i);
        double h = 1./(N+1);
        double hh = h*h;
        double start, stop;

        //M = midle diagonal
        //F = function evaluation
        //U = unknown function result

        double* M = new double[N+2];
        double* U = new double[N+2];
        double* F = new double[N+2];
        double* x = new double[N+2]; 

        x[1] = U[1] = F[1] = 0; x[N+1] = U[N+1] = F[N+1] = 0;
        M[1] = M[N] =  2.0;

        for (int i = 2; i< N+2; i++){
            x[i] = i*h;
            M[i] = (i + 1.0)/((double) i);
            F[i] = hh*funct(x[i]);
        }

        //start of algorithm
        start = clock();

        for(int i = 2; i < N+2 ; i++){
            F[i] = F[i] + F[i-1]/M[i-1] ;
        }

         U[N] = F[N]/M[N];

         for(int i = N-1 ; i > 1; i--){
             U[i] = (F[i] + U[i+1])/M[i];
         }
         stop = clock();
         //end of algorithm
         time_write << N <<"  "<<(double (stop - start) / CLOCKS_PER_SEC) << endl ;

         ofstream rel_file (rel_err_file.c_str());


         for(int i = 1; i < N +2; i++){
            double exact_val  = exact(x[i]);
            double Relerr = fabs( (exact_val - U[i]) / exact_val);

            rel_file << setw(15) << setprecision(8) << x[i];
            rel_file << setw(15) << setprecision(8) << U[i];
            rel_file << setw(15) << setprecision(8) << exact_val;
            rel_file << setw(15) << setprecision(8) << log10(Relerr) << endl;
            }

         rel_file.close();
         delete M; delete F; delete x; delete U;
    }
time_write.close();
return 0;
}


int main(int argc, char *argv[]){
    int exponent = atoi(argv[1]);
    int N = pow(10, exponent);

    double* U = new double[N+2];
    non_special_method(exponent);
    special_method(exponent);

    return 0;

}
