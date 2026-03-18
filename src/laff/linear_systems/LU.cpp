#include "laff/Laff.hpp"

namespace laff {
    bool lu_unb(Matrix& A) {
        if (A.m != A.n) return false;
        int n = A.m;
        for (int k = 0; k < n; k++) {
            double alpha_kk = A(k, k);
            if (alpha_kk == 0.0) return false;
            for (int i = k + 1; i < n; i++) A(i, k) /= alpha_kk;
            for (int j = k + 1; j < n; j++) {
                for (int i = k + 1; i < n; i++) {
                    A(i, j) -= A(i, k) * A(k, j);
                }
            }
        }
        return true;
    }

    bool solve_lu(Matrix& A, Matrix& b) {
        if (!trsv_ln(A, b)) return false;
        if (!trsv_un(A, b)) return false;
        return true;
    }
    // ... (rest of implementation)
}
