#include "laff/Laff.hpp"

namespace laff {
    bool ger(double alpha, const Matrix& x, const Matrix& y, Matrix& A) {
        if (x.m * x.n != A.m || y.m * y.n != A.n) return false;
        for (int j = 0; j < A.n; j++) {
            axpy(alpha * y.data[j], x, const_cast<Matrix&>(A).col(j));
        }
        return true;
    }

    bool syr_l(double alpha, const Matrix& x, Matrix& A) {
        if (x.m * x.n != A.m || A.m != A.n) return false;
        for (int j = 0; j < A.n; j++) {
            for (int i = j; i < A.m; i++) {
                A(i, j) += alpha * x.data[i] * x.data[j];
            }
        }
        return true;
    }
    // ... (rest of functions)
}
