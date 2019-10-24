#include "Matrix.h"
#include <vector>
#include <cmath>

using namespace mat_vec;
using namespace std;

Matrix::Matrix(size_t size, double value) : matSizeRow(size), matSizeCol(size), mat(new double * [size]) {
    for (size_t i = 0; i < matSizeRow; i++) {
        mat[i]  = new double[matSizeCol];
        for (size_t j = 0; j < matSizeCol; j++) {
            mat[i][j] = value;
        }
    }
}

Matrix Matrix::eye(size_t size) {
    Matrix res(size);
    for (size_t i = 0; i < res.matSizeRow; i++) {
        for (size_t j = 0; j < res.matSizeCol; j++) {
            if (i == j) {
                res.mat[i][j] = 1;
            }
        }
    }
    return res;
}

Matrix::Matrix(size_t rows, size_t cols, double value) : matSizeRow(rows), matSizeCol(cols), mat(new double * [rows]) {
    for (size_t i = 0; i < matSizeRow; i++) {
        mat[i] = new double[matSizeCol];
        for (size_t j = 0; j < matSizeCol; j++) {
            mat[i][j] = value;
        }
    }
}

Matrix::Matrix(const Matrix &src) : matSizeRow(src.matSizeRow), matSizeCol(src.matSizeCol), mat(new double*[matSizeRow]) {
    for (size_t i = 0; i < matSizeRow; i++) {
        mat[i] = new double[matSizeCol];
        for (size_t j = 0; j < matSizeCol; j++) {
            mat[i][j] = src.mat[i][j];
        }
    }
}

Matrix& Matrix::operator=(const Matrix &rhs) {
    for (size_t i = 0; i < matSizeRow; i++) {
        delete [] mat[i];
    }
    delete [] mat;
    matSizeRow = rhs.matSizeRow;
    matSizeCol = rhs.matSizeCol;
    mat = new double*[matSizeRow];
    for (size_t i = 0; i < matSizeRow; i++) {
        mat[i] = new double[matSizeCol];
        for (size_t j = 0; j < matSizeCol; j++) {
            mat[i][j] = rhs.mat[i][j];
        }
    }
    return *this;
}

std::pair<size_t, size_t> Matrix::shape() const {
    return std::pair<size_t, size_t>(matSizeRow,matSizeCol);
}

Matrix::~Matrix() {
    for (size_t i = 0; i < matSizeRow; i++) {
        delete [] mat[i];
    }
    delete [] mat;
    matSizeRow = 0;
    matSizeCol = 0;
}

void Matrix::reshape(size_t rows, size_t cols) {
    Matrix res(rows, cols, 0);

    size_t i = 0, i1 = 0, j = 0, j1 = 0;

    while ((i + 1) * (j + 1) <= matSizeRow * matSizeCol) {
        if (j == matSizeCol) {
            i++;
            j = 0;
        }
        if (j1 == cols) {
            i1++;
            j1 = 0;
        }
        res.get(i1, j1) = mat[i][j];
        j++;
        j1++;
    }

    *this = res;
}

double Matrix::get(size_t row, size_t col) const {
    return mat[row][col];
}

double& Matrix::get(size_t row, size_t col) {
    return mat[row][col];
}

Matrix Matrix::operator+(const Matrix &rhs) const {
    Matrix res(matSizeRow, matSizeCol);
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            res.mat[i][j] = mat[i][j] + rhs.mat[i][j];
        }
    }
    return res;
}

Matrix& Matrix::operator+=(const Matrix &rhs) {
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            mat[i][j] += rhs.mat[i][j];
        }
    }
    return *this;
}

Matrix Matrix::operator-(const Matrix &rhs) const {
    Matrix res(matSizeRow, matSizeCol);
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            res.mat[i][j] = mat[i][j] - rhs.mat[i][j];
        }
    }
    return res;
}

Matrix& Matrix::operator-=(const Matrix &rhs) {
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            mat[i][j] -= rhs.mat[i][j];
        }
    }
    return *this;
}

Matrix Matrix::operator*(const Matrix &rhs) const {
    Matrix res(matSizeRow, rhs.matSizeCol, 0);
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < rhs.matSizeCol; j++) {
            for (size_t n = 0; n < matSizeCol; n++) {
                res.mat[i][j] += mat[i][n] * rhs.mat[n][j];
            }
        }
    }
    return res;
}

Matrix& Matrix::operator*=(const Matrix &rhs) {
    Matrix res(matSizeRow, rhs.matSizeCol, 0);
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < rhs.matSizeCol; j++) {
            for (size_t n = 0; n < matSizeCol; n++) {
                res.mat[i][j] += mat[i][n] * rhs.mat[n][j];
            }
        }
    }
    *this = res;
    return *this;
}

Matrix Matrix::operator*(double k) const {
    Matrix res(matSizeRow, matSizeCol);
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            res.mat[i][j] = mat[i][j] * k;
        }
    }
    return res;
}

Matrix& Matrix::operator*=(double k) {
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            mat[i][j] *= k;
        }
    }
    return *this;
}

Matrix Matrix::operator/(double k) const {
    Matrix res(matSizeRow, matSizeCol);
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            res.mat[i][j] = mat[i][j] / k;
        }
    }
    return res;
}

Matrix& Matrix::operator/=(double k) {
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            mat[i][j] /= k;
        }
    }
    return *this;
}

void Matrix::transpose() {
    Matrix res(matSizeCol, matSizeRow);
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            res.mat[j][i] = mat[i][j];
        }
    }
    *this = res;
}

Matrix Matrix::transposed() const {
    Matrix res(matSizeCol, matSizeRow);
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            res.mat[j][i] = mat[i][j];
        }
    }
    return res;
}

bool Matrix::operator==(const mat_vec::Matrix &rhs) const {
    if (matSizeRow != rhs.matSizeRow || matSizeCol != rhs.matSizeCol) {
        return false;
    }
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            if (rhs.mat[i][j] != mat[i][j]) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::operator!=(const mat_vec::Matrix &rhs) const {
    return !(*this == rhs);
}

Vector Matrix::operator*(const mat_vec::Vector &vec) const {
    Vector res(matSizeRow);
    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            res[i] += mat[i][j] * vec[j];
        }
    }
    return res;
}

double calcDet(vector<vector<double>> &Matrix) {
    double det = 0;
    if (Matrix.size() == 1) {
        return Matrix[0][0];
    }
    else if (Matrix.size() == 2) {
        det = Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];
        return det;
    }
    else {
        for (int p = 0; p < Matrix[0].size(); p++) {
            vector<vector<double>> TempMatrix;
            for (int i = 1; i < Matrix.size(); i++) {
                vector<double> TempRow;
                for (int j = 0; j < Matrix[i].size(); j++) {
                    if (j != p) {
                        TempRow.push_back(Matrix[i][j]);
                    }
                }
                if (!TempRow.empty())
                    TempMatrix.push_back(TempRow);
            }
            det = det + Matrix[0][p] * pow(-1, p) * calcDet(TempMatrix);
        }
        return det;
    }
}

double Matrix::det() const {
    vector<vector<double>> matvec(matSizeRow, vector<double>(matSizeRow));
    for (int i = 0; i < matSizeRow; i++) {
        for (int j = 0; j < matSizeRow; ++j) {
            matvec[i][j] = get(i, j);
        }
    }
    return calcDet(matvec);
}

Matrix Matrix::inv() const {
    double d = this->det();
    Matrix res(matSizeRow);
    Matrix tmp(matSizeRow - 1);

    for (size_t i = 0; i < matSizeRow; i++) {
        for (size_t j = 0; j < matSizeCol; j++) {
            res.get(i, j) = (i + j) % 2 ? -1 : 1;
            size_t n = 0, m = 0, i1 = 0, j1 = 0;
            while ((i1 + 1) * (j1 + 1) <= (matSizeRow - 1) * (matSizeRow - 1)) {
                if (j1 == matSizeRow - 1) {
                    i1++;
                    j1 = 0;
                }
                if (m == matSizeRow) {
                    n++;
                    m = 0;
                }
                n += (n == i);
                m += (m == j);
                if (m == matSizeRow) {
                    n++;
                    m = 0;
                }
                n += (n == i);
                tmp.mat[i1][j1] = mat[n][m];
                m++;
                j1++;
            }
            res.get(i, j) *= tmp.det();
        }
    }

    res.transpose();
    res /= d;
    return res;
}

