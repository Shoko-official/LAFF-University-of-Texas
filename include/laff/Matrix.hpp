#ifndef MATRIX_HPP
#define MATRIX_HPP

/**
 * Basic header of matrices
 * Created the 2026/02/15 by Shoko_ofi
 * Edited the 2026/03/22 by Shoko (Added col/row methods)
 */

struct Matrix {
    int m;              // rows
    int n;              // cols
    int ldim;           // leading dimension
    double* data;
    bool owns_memory;

    Matrix(int rows = 0, int cols = 0, double val = 0.0);
    Matrix(double* ptr, int rows, int cols, int leading_dim);
    ~Matrix();

    Matrix(const Matrix& other);
    Matrix& operator=(const Matrix& other);

    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;

    // View methods for columns and rows
    Matrix col(int j);
    Matrix row(int i);

    Matrix slice(int row_start, int row_end, int col_start, int col_end);
};

#endif // MATRIX_HPP
