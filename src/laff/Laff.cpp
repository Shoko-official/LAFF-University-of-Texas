#include "laff/Laff.hpp"

namespace laff {
    bool fill(Matrix& A, double val) {
        for (int i = 0; i < A.m * A.n; i++) A.data[i] = val;
        return true;
    }

    bool zeros(Matrix& A) {
        return fill(A, 0.0);
    }

    bool scal_matrix(double alpha, Matrix& A) {
        for (int i = 0; i < A.m * A.n; i++) A.data[i] *= alpha;
        return true;
    }

    bool identity(Matrix& A) {
        if (A.m != A.n) return false;
        zeros(A);
        for (int i = 0; i < A.m; i++) A(i, i) = 1.0;
        return true;
    }

    bool diag(const Matrix& x, Matrix& A) {
        if (x.m * x.n != A.m || A.m != A.n) return false;
        zeros(A);
        for (int i = 0; i < A.m; i++) A(i, i) = x.data[i];
        return true;
    }

    bool transpose(const Matrix& A, Matrix& B) {
        if (A.m != B.n || A.n != B.m) return false;
        for (int i = 0; i < A.m; i++) {
            for (int j = 0; j < A.n; j++) {
                B(j, i) = A(i, j);
            }
        }
        return true;
    }

    bool lower_tri(Matrix& A) {
        for (int j = 0; j < A.n; j++) {
            for (int i = 0; i < j && i < A.m; i++) {
                A(i, j) = 0.0;
            }
        }
        return true;
    }

    bool upper_tri(Matrix& A) {
        for (int j = 0; j < A.n; j++) {
            for (int i = j + 1; i < A.m; i++) {
                A(i, j) = 0.0;
            }
        }
        return true;
    }

    bool symmetrize_from_lower(Matrix& A) {
        if (A.m != A.n) return false;
        for (int j = 0; j < A.n; j++) {
            for (int i = j + 1; i < A.m; i++) {
                A(j, i) = A(i, j);
            }
        }
        return true;
    }

    bool symmetrize_from_upper(Matrix& A) {
        if (A.m != A.n) return false;
        for (int j = 0; j < A.n; j++) {
            for (int i = j + 1; i < A.m; i++) {
                A(i, j) = A(j, i);
            }
        }
        return true;
    }

    bool scal_matrix(double alpha, Matrix& A) {
        for (int i = 0; i < A.m * A.n; i++) A.data[i] *= alpha;
        return true;
    }

    bool add_matrix(const Matrix& B, Matrix& A) {
        if (A.m != B.m || A.n != B.n) return false;
        for (int i = 0; i < A.m * A.n; i++) A.data[i] += B.data[i];
        return true;
    }
}
