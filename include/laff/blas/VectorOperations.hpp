#ifndef VECTOR_OPERATIONS_HPP
#define VECTOR_OPERATIONS_HPP

#include "../Matrix.hpp"

namespace laff {
    bool copy(const Matrix& x, Matrix& y);
    bool scal(double alpha, Matrix& x);
    bool axpy(double alpha, const Matrix& x, Matrix& y);
    bool dot(const Matrix& x, const Matrix& y, double& alpha);
    bool norm2(const Matrix& x, double& alpha);
}

#endif // VECTOR_OPERATIONS_HPP
