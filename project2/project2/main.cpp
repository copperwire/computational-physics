#include <iostream>
#include <fstream>
#include <armadillo>
#include <string>
#include <cmath>

using namespace std;

double V(double x){
    return x*x;
}

void max_elem(int i, int j, arma::mat A, int N){
    for(int k = 0; i < N+1 ;i++ ){
        for(int l = 0; j < N+1 ; j++ ){
            if(A[i, j] < A[k, l]){
                i = k;
                j = l;
            }

        }
    }
}

int main(int argc, char *argv[])
{
    double epsilon = 1e-8;
    double exponent = log10(2) ; // atoi(argv[1]);
    int rhomax = 4; //atoi(argv[2]);




    int N = pow(10, exponent);
    int num__iter = 0 ;
    int max_iter = 100;

    double h = rhomax /(N - 1);
    int i, j = 0;

    double nondiagconst = - 1/(h*h) ;
    double diagconst = - 2 * nondiagconst;

    arma::vec rho(N+1);
    arma::vec U(N+1);
    arma::mat A(N+1, N+1);

    //arma::mat S(N+1);

    for(int i = 0; i < N+1 ;i++ ){
        for(int j = 0; j < N+1 ; j++ ){
            if(i == j){
                A(i,j) = diagconst ;
                //cout <<A[i, j] <<endl;
            }
            else if(i == j-1  || i == j+1){
                    A(i, j) = nondiagconst;
                }
            else {};
        //end of j
        }

    //end of i
    }

    while (num__iter <= max_iter )
        // must add length condition too
    {

    }

    return 0;
}
