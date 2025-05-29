#include "gaussIntegral.h"

myMatrix GaussInt4Q4(const myMatrix& D, const double ksi, const double eta, const double weight, const double x[4], const double y[4]) {
    int J_1=1;
    myMatrix adJ(2,2);
    myMatrix B(3,8);
    myMatrix Bt(8,3);
    myMatrix G(2,4); 
    myMatrix numdaN(2,4);
    myMatrix K(8,8);
    //calculate Jacobian matrix
        adJ(1,1) = 0.25 * ((x[1] - x[0])*(1-eta) - (x[3] - x[2])*(1+eta));
        adJ(1,0) = -0.25 * ((x[3] - x[0])*(1-eta) + (x[2] - x[1])*(1+eta));
        adJ(0,1) = -0.25 * ((y[1] - y[0])*(1-ksi) - (y[3] - y[2])*(1+ksi));
        adJ(0,0) = 0.25 * ((y[3] - y[0])*(1-ksi) + (y[2] - y[1])*(1+ksi));
        J_1 = 1/(adJ(0,0) * adJ(1,1) - adJ(0,1) * adJ(1,0));
        G(0,0) = eta-1 ;
        G(0,1) = 1-eta ;
        G(0,2) = 1+eta ;
        G(0,3) = -eta-1 ;
        G(1,0) = ksi-1 ;
        G(1,1) = -ksi-1 ;
        G(1,2) = 1+ksi ;
        G(1,3) = 1-ksi ;
        numdaN=(adJ * G) * J_1; // N = J^-1 * G
        //calculate B matrix
        B(0,0) = numdaN(0,0);
        B(0,2) = numdaN(0,1);
        B(0,4) = numdaN(0,2);
        B(0,6) = numdaN(0,3);
        B(1,1) = numdaN(1,0);
        B(1,3) = numdaN(1,1);
        B(1,5) = numdaN(1,2);
        B(1,7) = numdaN(1,3);
        B(2,0) = numdaN(1,0);
        B(2,1) = numdaN(0,0);
        B(2,2) = numdaN(1,1);
        B(2,3) = numdaN(0,1);
        B(2,4) = numdaN(1,2);
        B(2,5) = numdaN(0,2);
        B(2,6) = numdaN(1,3);
        B(2,7) = numdaN(0,3);
        Bt = B.transpose();

        //calculate K matrix
        K = weight * (Bt * D * B) * J_1; // K = sum  B^T * D * B * J^-1

    return K;
}