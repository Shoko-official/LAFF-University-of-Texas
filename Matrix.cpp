#include "Matrix.hpp"

/*
 * Basic implementation of matrices
 * Created the 2026/02/15 by Shoko_ofi
 * Editted the 2026/03/08 by Shoko_ofi
 */

Matrix::Matrix(int rows, int cols, double val) : m(rows), n(cols), ldim(cols), owns_memory(true) {
    data = new double[m * n];
    for (int i = 0; i < m * n; ++i) {
        data[i] = val;
    }
}

Matrix::Matrix(double* ptr, int rows, int cols, int leading_dim) 
    : m(rows), n(cols), ldim(leading_dim), data(ptr), owns_memory(false) {}

Matrix::~Matrix() {
    if (owns_memory) {
        delete[] data;
    }
}

Matrix::Matrix(const Matrix& other) : m(other.m), n(other.n), ldim(other.n), owns_memory(true) {
    data = new double[m * n];
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            data[i * ldim + j] = other.data[i * other.ldim + j];
        }
    }
}

Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        if (owns_memory) {
            delete[] data;
        }
        m = other.m;
        n = other.n;
        ldim = other.n;
        owns_memory = true;
        data = new double[m * n];
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                data[i * ldim + j] = other.data[i * other.ldim + j];
            }
        }
    }
    return *this;
}

double& Matrix::operator()(int i, int j) {
    return data[i * ldim + j];
}

const double& Matrix::operator()(int i, int j) const {
    return data[i * ldim + j];
}

Matrix Matrix::slice(int row_start, int row_end, int col_start, int col_end) {
    return Matrix(data + (row_start * ldim) + col_start, row_end - row_start, col_end - col_start, ldim);
}
