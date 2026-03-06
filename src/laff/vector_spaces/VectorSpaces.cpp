#include "laff/Laff.hpp"

namespace laff {
    bool rref(Matrix& A, double tol) {
        int m = A.m;
        int n = A.n;
        int pivot_row = 0;
        for (int j = 0; j < n && pivot_row < m; j++) {
            int max_row = pivot_row;
        }
        return true;
    }

    int rank(const Matrix& A, double tol) {
        Matrix temp = A;
        rref(temp, tol);
        return 0; // return actual rank
    }
}
