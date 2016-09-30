#include <iostream>
#include <fstream>
#include <armadillo>
#include <string>
#include <cmath>

using namespace std;

double V(double x){
    return x*x;
}

double max_elem(int* i, int* j, arma::mat A, int N){
    // matrix is symmetric - needs only check the upper or lower part.
    double max;
    for(int k = 0; k < N +1; k++ ){
        for(int l = k+ 1; l < N +1; l++ ){

            double a_kl = fabs(A(k, l));

            if(a_kl > max){
                max = a_kl;
                *i = k;
                *j = l;
            }
        }
    }
    return max;
}

arma::mat sim_transform(int k, int l, arma::mat A, arma::mat R, int N){
    double s, c, t, tau;

    //A.print();

    if (A(k, l) != 0.0){
        tau = (A(l,l) - A(k,k))/(2.0*A(k,l));

        if(tau > 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
            //t = +tau - sqrt(1+tau*tau);
        }
        else{
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        c = 1/(sqrt(1.0 +t*t));
        s = c*t;
    }
    else{
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
    A(l,k) = 0.0;  // same here
    for ( int i = 0; i < N+1; i++ ) {
      if ( i != k && i != l ) {
        a_ik = A(i,k);
        a_il = A(i,l);
        A(i,k) = c*a_ik - s*a_il;
        A(k,i) = A(i,k);
        A(i,l) = c*a_il + s*a_ik;
        A(l,i) = A(i,l);
      }
      r_ik = R(i,k);
      r_il = R(i,l);

      R(i,k) = c*r_ik - s*r_il;
      R(i,l) = c*r_il + s*r_ik;
    }
    return A;
}

int main(int argc, char *argv[])
{
    double epsilon = 1e-8;
    double exponent = atoi(argv[1]);
    int rhomax = atoi(argv[2]);


    int N = 399; //pow(10, exponent);
    int num_iter = 0 ;
    int max_iter = 1e5;
    double max = 1;

    double h = rhomax /double (N);
    int i = 0;
    int j = 0;

    double nondiagconst = - 1.0/(h*h) ;
    double diagconst =  2.0/(h*h);

    arma::vec rho(N+1);
    arma::vec U(N+1);
    arma::mat A(N+1, N+1);
    arma::mat R(N+1, N+1); //eigenvec matrix
    arma::vec diag(N+1);
    arma::vec ref_eigval(N+1);

    A.zeros();
    R.eye();
    U.zeros();
    rho.zeros();

    for(int i = 0; i < N +1; i++){
        //cout << "i, h" << i << " " << h << endl;
        rho(i) = i*h;
    }

    for(int i = 0; i < N + 1 ;i++ ){
        for(int j = 0; j < N + 1 ; j++ ){
            if(i == j){
                A(i,j) = diagconst + V(rho[i]);
                //cout <<"element" <<A(i, j) << "index" << i << j << endl;
            }
            else if(i == j-1 || i == j+1){
                A(i, j) = nondiagconst ;
                //cout <<"element" << A[i, j] << "index" << i << j << endl;
                }
            else{
                A(i, j) = 0.0 ;
            }
        //end of j
        }

    //end of i
    }

    while (num_iter <= max_iter && fabs(max) > epsilon)
    {
        max = max_elem(&i, &j, A, N);
        A = sim_transform(i, j, A, R,N);
        num_iter++;
    }
    cout << "value of max matrix element | number of iterations" << "  " << max << "  " << num_iter << endl;

    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    diag = A.diag();

    //arma::eig_sym(ref_eigval, R, A);
    diag = arma::sort(diag);
    //ref_eigval = arma::sort(ref_eigval);
    //ref_eigval.print();
    //cout << "-------------------------------------------------------------------------------------------------------" << endl;
    diag.print();



    return 0;
}
