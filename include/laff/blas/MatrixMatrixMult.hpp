#ifndef MATRIX_MATRIX_MULT_HPP
#define MATRIX_MATRIX_MULT_HPP

#include "../Matrix.hpp"

namespace laff {
    bool ger(double alpha, const Matrix& x, const Matrix& y, Matrix& A);
    bool syr_l(double alpha, const Matrix& x, Matrix& A);
    bool gemm(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C);
    bool gemm_dot(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C);
    bool gemm_axpy(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C);
    bool gemm_row(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C);
    bool gemm_outer(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C);
    bool trsv_ln(const Matrix& L, Matrix& x);
    bool trsv_un(const Matrix& U, Matrix& x);
}

#endif // MATRIX_MATRIX_MULT_HPP
