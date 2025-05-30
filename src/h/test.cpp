#include <iostream>
#include <iomanip>
#include "MatrixMultiplication.h"

int main() {
    // 创建两个矩阵
    myMatrix A(2, 3);
    myMatrix B(3, 2);

    // 初始化A
    int val = 1;
    for (int i = 0; i < A.rowCount(); ++i)
        for (int j = 0; j < A.colCount(); ++j)
            A(i, j) = val++;

    // 初始化B
    val = 1;
    for (int i = 0; i < B.rowCount(); ++i)
        for (int j = 0; j < B.colCount(); ++j)
            B(i, j) = val++;

    // 输出A
    std::cout << "Matrix A:" << std::endl;
    for (int i = 0; i < A.rowCount(); ++i) {
        for (int j = 0; j < A.colCount(); ++j)
            std::cout << std::setw(5) << A(i, j) << " ";
        std::cout << std::endl;
    }

    // 输出B
    std::cout << "Matrix B:" << std::endl;
    for (int i = 0; i < B.rowCount(); ++i) {
        for (int j = 0; j < B.colCount(); ++j)
            std::cout << std::setw(5) << B(i, j) << " ";
        std::cout << std::endl;
    }

    // 矩阵乘法
    myMatrix C = A * B;
    std::cout << "A * B:" << std::endl;
    for (int i = 0; i < C.rowCount(); ++i) {
        for (int j = 0; j < C.colCount(); ++j)
            std::cout << std::setw(5) << C(i, j) << " ";
        std::cout << std::endl;
    }

    // 矩阵加法（A转置后与自身相加）
    myMatrix At = A.transpose();
    myMatrix D = At + At;
    std::cout << "A^T + A^T:" << std::endl;
    for (int i = 0; i < D.rowCount(); ++i) {
        for (int j = 0; j < D.colCount(); ++j)
            std::cout << std::setw(5) << D(i, j) << " ";
        std::cout << std::endl;
    }

    // 数乘
    myMatrix E = 2.0 * A;
    std::cout << "2 * A:" << std::endl;
    for (int i = 0; i < E.rowCount(); ++i) {
        for (int j = 0; j < E.colCount(); ++j)
            std::cout << std::setw(5) << E(i, j) << " ";
        std::cout << std::endl;
    }

    return 0;
}