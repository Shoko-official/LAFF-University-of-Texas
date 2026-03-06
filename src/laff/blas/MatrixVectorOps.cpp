#include "laff/Laff.hpp"

namespace laff {
    bool gemv_dot(const Matrix& A, const Matrix& x, Matrix& y) {
        if (x.m * x.n != A.n || y.m * y.n != A.m) return false;
        for (int i = 0; i < A.m; i++) {
            Matrix row = const_cast<Matrix&>(A).row(i);
            dot(row, x, y.data[i]);
        }
        return true;
    }

    bool gemv_axpy(const Matrix& A, const Matrix& x, Matrix& y) {
        if (x.m * x.n != A.n || y.m * y.n != A.m) return false;
        zeros(y);
        for (int j = 0; j < A.n; j++) {
            Matrix col = const_cast<Matrix&>(A).col(j);
            axpy(x.data[j], col, y);
        }
        return true;
    }
}
