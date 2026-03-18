#ifndef CHOLESKY_HPP
#define CHOLESKY_HPP

#include "../Matrix.hpp"

namespace laff {
    bool chol(Matrix& A);
    bool solve_chol(const Matrix& L, Matrix& b);
}

#endif // CHOLESKY_HPP
