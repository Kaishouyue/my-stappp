#include "gaussIntegral.h"
#include <iostream>
#include <iomanip>

myMatrix GaussInt4Q4(const myMatrix& D, const double ksi, const double eta, const double weight, const double x[4], const double y[4]) {
    double J_1=1.0;
    myMatrix adJ(2,2);
    myMatrix B(3,8);
    myMatrix Bt(8,3);
    myMatrix G(2,4); 
    myMatrix numdaN(2,4);
    myMatrix K(8,8);
    myMatrix J(2,2);
    myMatrix XY(4,2);
    //calculate XY matrix
    for(int i = 0; i < 4; ++i) {
        XY(i, 0) = x[i];
        XY(i, 1) = y[i];
    }

    //calculate Jacobian matrix
    //    adJ(1,1) = 0.25 * ((x[1] - x[0])*(1-eta) - (x[3] - x[2])*(1+eta));
    //    adJ(1,0) = -0.25 * ((x[3] - x[0])*(1-eta) + (x[2] - x[1])*(1+eta));
    //   adJ(0,1) = -0.25 * ((y[1] - y[0])*(1-ksi) - (y[3] - y[2])*(1+ksi));
    //   adJ(0,0) = 0.25 * ((y[3] - y[0])*(1-ksi) + (y[2] - y[1])*(1+ksi));
    
        G(0,0) = eta-1 ;
        G(0,1) = 1-eta ;
        G(0,2) = 1+eta ;
        G(0,3) = -eta-1 ;
        G(1,0) = ksi-1 ;
        G(1,1) = -ksi-1 ;
        G(1,2) = 1+ksi ;
        G(1,3) = 1-ksi ;
        G = G * 0.25;
        J= G * XY; // J = G * XY
    //calculate adJ matrix
    adJ(0,0) = J(1,1);
    adJ(0,1) = -J(0,1);
    adJ(1,0) = -J(1,0);
    adJ(1,1) = J(0,0);
    //calculate J_1
    J_1 = 1/(J(0,0) * J(1,1) - J(0,1) * J(1,0)); // J_1 = 1/|J|
    numdaN= adJ * G;
    numdaN= J_1 * numdaN; // N = J^-1 * G
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
#ifdef _DEBUG_
    // 输出J矩阵
    std::cout << "adJ matrix:" << std::endl;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << std::setw(12) << J(i, j) << " ";
        }
        std::cout << std::endl;
    }
    // 输出G矩阵
    std::cout << "G matrix:" << std::endl;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::cout << std::setw(12) << G(i, j) << " ";
        }
        std::cout << std::endl;
    }
    // 输出numdaN矩阵
    std::cout << "numdaN matrix:" << std::endl;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::cout << std::setw(12) << numdaN(i, j) << " ";
        }
        std::cout << std::endl;
    }
    // 输出B矩阵
    std::cout << "B matrix:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 8; ++j) {
            std::cout << std::setw(12) << B(i, j) << " ";
        }
        std::cout << std::endl;
    }
#endif
        //calculate K matrix
        K = weight * (Bt * D * B) * (1/J_1); // K = sum w B^T * D * B * |J|
       //  K = K;
    return K;
}

void GaussStress4Q4_Center(const myMatrix& D, const double x[4], const double y[4], 
                           const myMatrix& u, double* stress)
{
    // 单点高斯积分，中心点
    double ksi = 0.0, eta = 0.0;
    myMatrix adJ(2,2);
    myMatrix B(3,8);
    myMatrix G(2,4); 
    myMatrix numdaN(2,4);
    myMatrix J(2,2);
    myMatrix XY(4,2);

    // 计算XY矩阵
    for(int i = 0; i < 4; ++i) {
        XY(i, 0) = x[i];
        XY(i, 1) = y[i];
    }

    // 计算G矩阵
    G(0,0) = eta-1 ;
    G(0,1) = 1-eta ;
    G(0,2) = 1+eta ;
    G(0,3) = -eta-1 ;
    G(1,0) = ksi-1 ;
    G(1,1) = -ksi-1 ;
    G(1,2) = 1+ksi ;
    G(1,3) = 1-ksi ;
    G = G * 0.25;
    J = G * XY;

    // 计算adJ
    adJ(0,0) = J(1,1);
    adJ(0,1) = -J(0,1);
    adJ(1,0) = -J(1,0);
    adJ(1,1) = J(0,0);

    // 计算Jacobian行列式
    double detJ = J(0,0) * J(1,1) - J(0,1) * J(1,0);
    double J_1 = 1.0 / detJ;

    // 计算numdaN
    numdaN = adJ * G;
    numdaN = J_1 * numdaN;

    // 组装B矩阵
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

    // 直接用u参与运算
    myMatrix strain = B * u; // 3x1
    myMatrix sigma = D * strain; // 3x1

    // 输出到stress数组
    for(int i=0; i<3; ++i) stress[i] = sigma(i,0);
}