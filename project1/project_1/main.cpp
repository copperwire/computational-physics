#include <iostream>
#include <math.h>
#include <armadillo>
#include <stdlib.h>

using namespace std;

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    double h = 1./(N+1);

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


    /* Remember to write to file
     * Implement the exact function in python - faster
     * NB!
     * The function needs to be EXACTLY 0 at boundary conditions DIRICHLET MOTHERFUCKER
    c*/

    for(int i = 0; i<N+2; i++){
        A[i] = -1.0;
        B[i] = -1.0;
        M[i] = 2.0;
        F[i] = h*h*100*exp(-10* i*h);

    }

    for(int i = 2; i < N +2; i++){
        M[i] = M[i] - (A[i] * B[i-1])/M[i-1];
        F[i] = F[i] - (A[i] * F[i-1])/M[i-1];
    }

    for(int i = N; i > 1; i--){
        U[i - 1] = (F[i-1] -  C[i-1] * U[i])/B[i-1] ;
    }



}
