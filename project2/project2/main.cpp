#include <iostream>
#include <fstream>
#include <armadillo>
#include <string>
#include <cmath>
#include <stdlib.h>

using namespace std;

int test(arma::mat R, arma::mat A, arma::mat B, int N){

    arma::vec v(N+1);
    arma::vec w(N+1);
    arma::vec diag(N+1);
    arma::vec ref_eigval(N+1);


    int rand_elem = rand() % N ;
    int other_elem = rand() % N ;
    v = R(rand_elem);
    w = R(other_elem);

    if(dot(v, w) > 1e-6 && rand_elem != other_elem){
        cout<< "Orthogonality not preserved" << endl;
        cout <<"dot prod" << "   " << dot(v, w) << "  " <<rand_elem << other_elem << endl;
    }
    else if( fabs((dot(v, w) -1)) > 1e-6 && rand_elem == other_elem){
        cout<< "Orthogonality not preserved" << endl;
        cout << "dot prod" << "   " << dot(v, w) << "  " <<rand_elem << other_elem << endl;
    }

    arma::eig_sym(ref_eigval, R, B);
    diag = A.diag();
    diag = arma::sort(diag);
    ref_eigval = arma::sort(ref_eigval);

    for(int i = 0; i<3 ; i++){
        if(fabs(diag(i) - ref_eigval(i)) > 1e-4){
            cout << "coherence with eig_sym breached" << endl;
            cout << "value of difference:    " <<fabs(diag(i) - ref_eigval(i)) << endl;
        }
    }
}

double V(void*, double omega, double x){
    return omega*omega*x*x;
}
double other_V(void*, double omega, double x){
    return omega*omega*x*x  + 1/x ;
}

double max_elem(int* i, int* j, arma::mat A, int N){
    // matrix is symmetric - needs only check the upper or lower part.
    double max;
    for(int k = 0; k < N +1; k++ ){
        for(int l = k+1; l < N +1; l++ ){
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

int do_method(int N, int rhomax, double (*potential)(void*, double, double), double epsilon, void* context){

    double omega_vals [4] = {0.01, 0.5, 1, 5};
    for(int y = 0; y<4 ; y++){

        int num_iter = 0 ;
        int max_iter = N*N*N;
        double max = 1;

        double h = rhomax /double (N+1);
        int i = 0;
        int j = 0;

        double nondiagconst = - 1.0/(h*h) ;
        double diagconst =  2.0/(h*h);

        arma::vec rho(N+1);
        arma::vec U(N+1);
        arma::mat A(N+1, N+1);
        arma::mat R(N+1, N+1); //eigenvec matrix
        arma::vec diag(N+1);

        A.zeros();
        R.eye();
        U.zeros();
        rho.zeros();
        diag.zeros();

        for(int i = 0; i < N +1; i++){
            //cout << "i, h" << i << " " << h << endl;
            rho(i) = (i+1)*h;
        }

        double omega = omega_vals[y];
        for(int i = 0; i < N + 1 ;i++ ){
            for(int j = 0; j < N + 1 ; j++ ){
                if(i == j){
                    A(i,j) = diagconst + potential(context, omega, rho(i));
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

        arma::mat B(N+1, N+1);
        B = A;

        while (num_iter <= max_iter && fabs(max) > epsilon)
        {
            max = max_elem(&i, &j, A, N);
            A = sim_transform(i, j, A, R,N);
            num_iter++;
        }

        cout << "value of max matrix element | number of iterations | value of OMEGA" << "  " << max << "  " << num_iter << "  " << omega <<  endl;

        cout << "-------------------------------------------------------------------------------------------------------" << endl;
        diag = A.diag();
        diag = arma::sort(diag);

        test(R, A, B, N);

        for(int i = 0; i<3 ; i++){
            cout << diag(i) << endl;
        }

        /*  WRITING TO FILE GOES HERE
         * R = matrix of eigenvectors NB UNSORTED (eigenstates of the hamiltonian)
         * diag = Eigenvalues - sorted after line 177
         * */
        cout << "-------------------------------------------------------------------------------------------------------" << endl;
        }
     return 0;
}

int main(int argc, char *argv[])
{
    double epsilon = 1e-8;
    double exponent = atoi(argv[1]);
    int rhomax = atoi(argv[2]);
    int N = 20; //pow(10, exponent);

    // needs to expand so do_method takes filename as argument
    do_method(N, rhomax, &V, epsilon, 0);
    do_method(N, rhomax, &other_V, epsilon, 0);

}

