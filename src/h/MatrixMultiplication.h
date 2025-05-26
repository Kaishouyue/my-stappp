#ifndef MATRIX_MULTIPLICATION_H
#define MATRIX_MULTIPLICATION_H

#include <vector>
#include <stdexcept>

class myMatrix {
private:
    std::vector<std::vector<double>> data;
    int rows, cols;

public:
    myMatrix(int r, int c) : rows(r), cols(c), data(r, std::vector<double>(c, 0.0)) {}

    int rowCount() const { return rows; }
    int colCount() const { return cols; }

    double& operator()(int r, int c) {
        if (r < 0 || r >= rows || c < 0 || c >= cols)
            throw std::out_of_range("Matrix index out of range");
        return data[r][c];
    }

    double operator()(int r, int c) const {
        if (r < 0 || r >= rows || c < 0 || c >= cols)
            throw std::out_of_range("Matrix index out of range");
        return data[r][c];
    }

    myMatrix operator*(const myMatrix& other) const {
        if (cols != other.rows)
            throw std::invalid_argument("矩阵维度不匹配，无法相乘。");

        myMatrix result(rows, other.cols);
        for (int i = 0; i < rows; ++i)
            for (int k = 0; k < cols; ++k)
                for (int j = 0; j < other.cols; ++j)
                    result(i,j) += data[i][k] * other.data[k][j];

        return result;
    }

    myMatrix operator+(const myMatrix& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("矩阵维度不匹配，无法相加。");

        myMatrix result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result(i,j) = data[i][j] + other.data[i][j];

        return result;
    }

    myMatrix transpose() const {
        myMatrix result(cols, rows);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result(j,i) = data[i][j];
        return result;
    }

    myMatrix operator*(double scalar) const {
        myMatrix result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result(i,j) = data[i][j] * scalar;
        return result;
    }

    friend myMatrix operator*(double scalar, const myMatrix& mat) {
        return mat * scalar;  
    }
};

#endif // MATRIX_MULTIPLICATION_H
