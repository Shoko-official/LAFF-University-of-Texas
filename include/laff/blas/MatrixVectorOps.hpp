#ifndef MATRIX_VECTOR_OPS_HPP
#define MATRIX_VECTOR_OPS_HPP

#include "../Matrix.hpp"

namespace laff {
    bool gemv_dot(const Matrix& A, const Matrix& x, Matrix& y);
    bool gemv_axpy(const Matrix& A, const Matrix& x, Matrix& y);
    bool gemv_t_dot(const Matrix& A, const Matrix& x, Matrix& y);
    bool gemv_t_axpy(const Matrix& A, const Matrix& x, Matrix& y);
    bool trmv_ln(const Matrix& L, Matrix& x);
    bool trmv_un(const Matrix& U, Matrix& x);
    bool symv_l(const Matrix& A, const Matrix& x, Matrix& y);
}

#endif // MATRIX_VECTOR_OPS_HPP
