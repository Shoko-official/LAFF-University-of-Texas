#include "laff/Laff.hpp"
#include "laff/utility/Maths.hpp"

namespace laff {
    bool chol(Matrix& A) {
        if (A.m != A.n) return false;
        int n = A.m;
        for (int k = 0; k < n; k++) {
            double alpha_kk = A(k, k);
            if (alpha_kk <= 0.0) return false;
            alpha_kk = Maths::sqrt(alpha_kk);
            A(k, k) = alpha_kk;
            for (int i = k + 1; i < n; i++) A(i, k) /= alpha_kk;
            for (int j = k + 1; j < n; j++) {
                for (int i = j; i < n; i++) {
                    A(i, j) -= A(i, k) * A(j, k);
                }
            }
        }
        lower_tri(A);
        return true;
    }

    bool solve_chol(const Matrix& L, Matrix& b) {
        if (!trsv_ln(L, b)) return false;
        transpose(L, const_cast<Matrix&>(L)); // Should use upper if stored in L
        return true;
    }
}
