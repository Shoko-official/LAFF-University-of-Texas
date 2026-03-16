#include "Laff.hpp"
#include "Utility/Maths.hpp"

/*
* Implemented by Shoko on 2026-03-06
*/

namespace laff {
    bool copy(const Matrix& x, Matrix& y) {
        if ((x.m != 1 && x.n != 1) || (y.m != 1 && y.n != 1)) {
            return false;
        }

        if (x.m * x.n != y.m * y.n) {
            return false;
        }

        int size = x.m * x.n;
        for (int i = 0; i < size; i++) {
            double val = (x.m == 1) ? x(0, i) : x(i, 0);
            if (y.m == 1) y(0, i) = val;
            else y(i, 0) = val;
        }
           return true;
    }


    bool scal(double alpha, Matrix& x){
        if (x.m != 1 && x.n != 1) {
            return false;
        }

        int size = x.m * x.n;
        for (int i = 0; i < size; i++) {
            if (x.m == 1) x(0, i) *= alpha;
            else x(i, 0) *= alpha;
        }

        return true;
    }

    bool axpy(double alpha, const Matrix& x, Matrix& y) {
        if ((x.m != 1 && x.n != 1) || (y.m != 1 && y.n != 1)) {
            return false;
        }
        if (x.m * x.n != y.m * y.n) {
            return false;
        }
        int size = x.m * x.n;

        for (int i = 0; i < size; i++) {
            double x_val = (x.m == 1) ? x(0, i) : x(i, 0);
            if (y.m == 1) y(0, i) = alpha * x_val + y(0, i);
            else y(i, 0) = alpha * x_val + y(i, 0);
        }
        return true;
    }

    bool dot(const Matrix& x, const Matrix& y, double& alpha){
        if ((x.m != 1 && x.n != 1) || (y.m != 1 && y.n != 1)) {
            return false;
        }
        if (x.m * x.n != y.m * y.n) {
            return false;
        }
        int size = x.m * x.n;
        alpha = 0.0;
        
        for (int i = 0; i < size; i++){
            double x_val = (x.m == 1) ? x(0, i) : x(i, 0);
            double y_val = (y.m == 1) ? y(0, i) : y(i, 0);
            alpha += x_val * y_val;
        }
        return true;
    }

    bool norm2(const Matrix& x, double& alpha) {
        if (x.m != 1 && x.n != 1) {
            return false;
        }

        if (!dot(x, x, alpha)) {
            return false;
        }

        alpha = Maths::sqrt(alpha);
        return true;
    }

    bool zeros(Matrix& A) {
        for (int i = 0; i < A.m; i++) {
            for (int j = 0; j < A.n; j++) {
                A(i, j) = 0.0;
            }
        }
        return true;
    }

    bool identity(Matrix& A) {
        if (A.m != A.n) {
            return false;
        }

        zeros(A);

        /* Set the diagonal elements to 1.0 */
        for (int i = 0; i < A.m; i++) {
            A(i, i) = 1.0;
        }
        return true;
    }

    bool diag(const Matrix& x, Matrix& A) {
        if (x.m != 1 && x.n != 1) {
            return false;
        }
        int min_dim = (A.m < A.n) ? A.m : A.n;
        if (x.m * x.n < min_dim) {
            return false;
        }

        zeros(A);
        for (int i = 0; i < min_dim; i++) {
            double val = (x.m == 1) ? x(0, i) : x(i, 0);
            A(i, i) = val;
        }
        return true;
    }

    bool transpose(const Matrix& A, Matrix& B) {
        if (A.m != B.n || A.n != B.m) {
            return false;
        }

        for (int i = 0; i < A.m; i++) {
            for (int j = 0; j < A.n; j++) {
                B(j, i) = A(i, j);
            }
        }
        return true;
    }

    bool lower_tri(Matrix& A) {
        for (int i = 0; i < A.m; i++) {
            for (int j = i + 1; j < A.n; j++) {
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
            for (int i = 0; i < j; i++) {
                A(i, j) = A(j, i);
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
        for (int i = 0; i < A.m; i++) {
            for (int j = 0; j < A.n; j++) {
                A(i, j) *= alpha;
            }
        }
        return true;
    }

    bool add_matrix(const Matrix& B, Matrix& A) {
        if (A.m != B.m || A.n != B.n) return false;
        for (int i = 0; i < A.m; i++) {
            for (int j = 0; j < A.n; j++) {
                A(i, j) += B(i, j);
            }
        }
        return true;
    }

    bool gemv_dot(const Matrix& A, const Matrix& x, Matrix& y) {
        // Validation: A (m x n), x (n x 1 or 1 x n), y (m x 1 or 1 x m)
        if (x.m * x.n != A.n || y.m * y.n != A.m) return false;

        for (int i = 0; i < A.m; i++) {
            Matrix row = const_cast<Matrix&>(A).slice(i, i + 1, 0, A.n);
            double rho;
            dot(row, x, rho);
            
            // Handle y as either row or column vector
            if (y.m == 1) y(0, i) += rho;
            else y(i, 0) += rho;
        }
        return true;
    }

    bool gemv_axpy(const Matrix& A, const Matrix& x, Matrix& y) {
        // Validation: A (m x n), x (n x 1 or 1 x n), y (m x 1 or 1 x m)
        if (x.m * x.n != A.n || y.m * y.n != A.m) return false;

        for (int j = 0; j < A.n; j++) {
            Matrix col = const_cast<Matrix&>(A).slice(0, A.m, j, j + 1);
            double alpha = (x.m == 1) ? x(0, j) : x(j, 0);
            axpy(alpha, col, y);
        }
        return true;
    }

    bool gemv_t_dot(const Matrix& A, const Matrix& x, Matrix& y) {
        // Validation: A (m x n), x (m x 1 or 1 x m), y (n x 1 or 1 x n)
        if (x.m * x.n != A.m || y.m * y.n != A.n) return false;

        for (int j = 0; j < A.n; j++) {
            Matrix col = const_cast<Matrix&>(A).slice(0, A.m, j, j + 1);
            double rho;
            dot(col, x, rho);
            
            if (y.m == 1) y(0, j) += rho;
            else y(j, 0) += rho;
        }
        return true;
    }

    bool gemv_t_axpy(const Matrix& A, const Matrix& x, Matrix& y) {
        // Validation: A (m x n), x (m x 1 or 1 x m), y (n x 1 or 1 x n)
        if (x.m * x.n != A.m || y.m * y.n != A.n) return false;

        for (int i = 0; i < A.m; i++) {
            Matrix row = const_cast<Matrix&>(A).slice(i, i + 1, 0, A.n);
            double alpha = (x.m == 1) ? x(0, i) : x(i, 0);
            axpy(alpha, row, y);
        }
        return true;
    }

    bool trmv_ln(const Matrix& L, Matrix& x) {
        // Validation: L (n x n), x (n x 1 or 1 x n)
        if (L.m != L.n || x.m * x.n != L.m) return false;

        // Use a temporary to avoid overwriting elements before they are used
        Matrix x_copy(x.m, x.n);
        copy(x, x_copy);
        zeros(x);

        for (int i = 0; i < L.m; i++) {
            for (int j = 0; j <= i; j++) {
                double val = (x_copy.m == 1) ? x_copy(0, j) : x_copy(j, 0);
                if (x.m == 1) x(0, i) += L(i, j) * val;
                else x(i, 0) += L(i, j) * val;
            }
        }
        return true;
    }

    bool trmv_un(const Matrix& U, Matrix& x) {
        // Validation: U (n x n), x (n x 1 or 1 x n)
        if (U.m != U.n || x.m * x.n != U.m) return false;

        Matrix x_copy(x.m, x.n);
        copy(x, x_copy);
        zeros(x);

        for (int i = 0; i < U.m; i++) {
            for (int j = i; j < U.n; j++) {
                double val = (x_copy.m == 1) ? x_copy(0, j) : x_copy(j, 0);
                if (x.m == 1) x(0, i) += U(i, j) * val;
                else x(i, 0) += U(i, j) * val;
            }
        }
        return true;
    }

    bool symv_l(const Matrix& A, const Matrix& x, Matrix& y) {
        // Validation: A (n x n), x (n x 1 or 1 x n), y (n x 1 or 1 x n)
        if (A.m != A.n || x.m * x.n != A.m || y.m * y.n != A.m) return false;

        for (int i = 0; i < A.m; i++) {
            double sum = 0.0;
            for (int j = 0; j < A.n; j++) {
                double x_val = (x.m == 1) ? x(0, j) : x(j, 0);
                if (j <= i) {
                    sum += A(i, j) * x_val;
                } else {
                    sum += A(j, i) * x_val; // Symmetric part
                }
            }
            if (y.m == 1) y(0, i) += sum;
            else y(i, 0) += sum;
        }
        return true;
    }

    bool ger(double alpha, const Matrix& x, const Matrix& y, Matrix& A) {
        if (x.m * x.n != A.m || y.m * y.n != A.n) return false;

        for (int j = 0; j < A.n; j++) {
            double y_val = (y.m == 1) ? y(0, j) : y(j, 0);
            for (int i = 0; i < A.m; i++) {
                double x_val = (x.m == 1) ? x(0, i) : x(i, 0);
                A(i, j) += alpha * x_val * y_val;
            }
        }
        return true;
    }

    bool syr_l(double alpha, const Matrix& x, Matrix& A) {
        if (A.m != A.n || x.m * x.n != A.m) return false;

        for (int j = 0; j < A.n; j++) {
            double x_j = (x.m == 1) ? x(0, j) : x(j, 0);
            for (int i = j; i < A.m; i++) {
                double x_i = (x.m == 1) ? x(0, i) : x(i, 0);
                A(i, j) += alpha * x_i * x_j;
            }
        }
        return true;
    }

    bool gemm(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C) {
        return gemm_axpy(alpha, A, B, beta, C);
    }

    bool gemm_dot(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C) {
        if (A.m != C.m || B.n != C.n || A.n != B.m) return false;

        for (int i = 0; i < C.m; i++) {
            for (int j = 0; j < C.n; j++) {
                Matrix row_A = const_cast<Matrix&>(A).slice(i, i + 1, 0, A.n);
                Matrix col_B = const_cast<Matrix&>(B).slice(0, B.m, j, j + 1);
                double rho;
                dot(row_A, col_B, rho);
                C(i, j) = alpha * rho + beta * C(i, j);
            }
        }
        return true;
    }

    bool gemm_axpy(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C) {
        if (A.m != C.m || B.n != C.n || A.n != B.m) return false;

        scal_matrix(beta, C);

        // Compute C := alpha * A * B + C
        for (int j = 0; j < C.n; j++) {
            for (int k = 0; k < A.n; k++) {
                double b_val = B(k, j);
                double scaled_alpha = alpha * b_val;
                for (int i = 0; i < A.m; i++) {
                    C(i, j) += scaled_alpha * A(i, k);
                }
            }
        }
        return true;
    }

    bool gemm_row(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C) {
        if (A.m != C.m || B.n != C.n || A.n != B.m) return false;

        scal_matrix(beta, C);

        for (int i = 0; i < C.m; i++) {
            for (int k = 0; k < A.n; k++) {
                double a_val = A(i, k);
                double scaled_alpha = alpha * a_val;
                for (int j = 0; j < C.n; j++) {
                    C(i, j) += scaled_alpha * B(k, j);
                }
            }
        }
        return true;
    }

    bool gemm_outer(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C) {
        if (A.m != C.m || B.n != C.n || A.n != B.m) return false;

        scal_matrix(beta, C);

        for (int k = 0; k < A.n; k++) {
            Matrix col_A = const_cast<Matrix&>(A).slice(0, A.m, k, k + 1);
            Matrix row_B = const_cast<Matrix&>(B).slice(k, k + 1, 0, B.n);
            ger(alpha, col_A, row_B, C);
        }
        return true;
    }
}