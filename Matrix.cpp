#include "Matrix.hpp"

/*
 * Basic implementation of matrices
 * Corrected by Shoko on 2026-03-06
 */

Matrix::Matrix(int rows, int cols) : m(rows), n(cols), data(rows * cols) {}

double &Matrix::operator()(int i, int j) {
    return data[i * n + j];
}

const double& Matrix::operator()(int i, int j) const {
    return data[i * n + j];
}