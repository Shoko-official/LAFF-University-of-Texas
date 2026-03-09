#ifndef VECTOR_SPACES_HPP
#define VECTOR_SPACES_HPP

#include "../Matrix.hpp"

namespace laff {
    bool rref(Matrix& A, double tol = 1e-9);
    int rank(const Matrix& A, double tol = 1e-9);
    bool is_linearly_independent(const Matrix& A, double tol = 1e-9);
    bool is_spanning(const Matrix& A, double tol = 1e-9);
    bool is_basis(const Matrix& A, double tol = 1e-9);
}

#endif // VECTOR_SPACES_HPP
