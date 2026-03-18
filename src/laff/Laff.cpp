#include "Laff.hpp"
#include "Utility/Maths.hpp"
#include <iostream>
#include <iomanip>

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

    bool trsv_ln(const Matrix& L, Matrix& x) {
        if (L.m != L.n || x.m * x.n != L.m) return false;

        int n = L.m;
        for (int i = 0; i < n; i++) {
            double& x_i = (x.m == 1) ? x(0, i) : x(i, 0);
            for (int j = 0; j < i; j++) {
                double x_j = (x.m == 1) ? x(0, j) : x(j, 0);
                x_i -= L(i, j) * x_j;
            }
        }
        return true;
    }

    bool trsv_un(const Matrix& U, Matrix& x) {
        if (U.m != U.n || x.m * x.n != U.m) return false;

        int n = U.m;
        for (int i = n - 1; i >= 0; i--) {
            double& x_i = (x.m == 1) ? x(0, i) : x(i, 0);
            for (int j = i + 1; j < n; j++) {
                double x_j = (x.m == 1) ? x(0, j) : x(j, 0);
                x_i -= U(i, j) * x_j;
            }
            if (U(i, i) == 0.0) return false;
            x_i /= U(i, i);
        }
        return true;
    }

    bool lu_unb(Matrix& A) {
        if (A.m != A.n) return false;

        int n = A.m;
        for (int k = 0; k < n; k++) {
            double alpha_kk = A(k, k);
            if (alpha_kk == 0.0) return false;

            for (int i = k + 1; i < n; i++) {
                A(i, k) /= alpha_kk;
            }

            for (int j = k + 1; j < n; j++) {
                for (int i = k + 1; i < n; i++) {
                    A(i, j) -= A(i, k) * A(k, j);
                }
            }
        }
        return true;
    }

    bool solve_lu(Matrix& A, Matrix& b) {
        if (!lu_unb(A)) return false;
        if (!trsv_ln(A, b)) return false;
        if (!trsv_un(A, b)) return false;
        return true;
    }

    bool iamax(const Matrix& x, int& index) {
        if (x.m != 1 && x.n != 1) return false;
        int n = x.m * x.n;
        index = 0;
        double max_val = -1.0;
        for (int i = 0; i < n; i++) {
            double val = (x.m == 1) ? x(0, i) : x(i, 0);
            double abs_val = (val < 0) ? -val : val;
            if (abs_val > max_val) {
                max_val = abs_val;
                index = i;
            }
        }
        return true;
    }

    bool swap_rows(Matrix& A, int i, int j) {
        if (i == j) return true;
        if (i < 0 || i >= A.m || j < 0 || j >= A.m) return false;
        for (int k = 0; k < A.n; k++) {
            double tmp = A(i, k);
            A(i, k) = A(j, k);
            A(j, k) = tmp;
        }
        return true;
    }

    bool lu_piv(Matrix& A, Matrix& p) {
        if (A.m != A.n || p.m * p.n != A.m) return false;
        int n = A.m;
        // Initialize permutation vector
        for (int i = 0; i < n; i++) {
            if (p.m == 1) p(0, i) = (double)i;
            else p(i, 0) = (double)i;
        }

        for (int k = 0; k < n; k++) {
            // Find pivot in column k, from row k to n-1
            Matrix col_slice = A.slice(k, n, k, k + 1);
            int pivot_idx_local;
            iamax(col_slice, pivot_idx_local);
            int pivot_idx = k + pivot_idx_local;

            // Record pivot
            if (p.m == 1) p(0, k) = (double)pivot_idx;
            else p(k, 0) = (double)pivot_idx;

            // Swap rows in A
            swap_rows(A, k, pivot_idx);

            double alpha_kk = A(k, k);
            if (alpha_kk == 0.0) return false;

            // Scale column
            for (int i = k + 1; i < n; i++) {
                A(i, k) /= alpha_kk;
            }

            // Rank-1 update
            for (int j = k + 1; j < n; j++) {
                for (int i = k + 1; i < n; i++) {
                    A(i, j) -= A(i, k) * A(k, j);
                }
            }
        }
        return true;
    }

    bool apply_piv(const Matrix& p, Matrix& b) {
        if (p.m * p.n != b.m * b.n) return false;
        int n = b.m * b.n;
        for (int k = 0; k < n; k++) {
            int pivot_idx = (int)((p.m == 1) ? p(0, k) : p(k, 0));
            // Swap elements in b
            if (b.m == 1) {
                double tmp = b(0, k);
                b(0, k) = b(0, pivot_idx);
                b(0, pivot_idx) = tmp;
            } else {
                double tmp = b(k, 0);
                b(k, 0) = b(pivot_idx, 0);
                b(pivot_idx, 0) = tmp;
            }
        }
        return true;
    }

    bool solve_lu_piv(Matrix& A, const Matrix& p, Matrix& b) {
        if (!apply_piv(p, b)) return false;
        if (!trsv_ln(A, b)) return false;
        if (!trsv_un(A, b)) return false;
        return true;
    }

    bool inv(const Matrix& A, Matrix& Ainv) {
        if (A.m != A.n || Ainv.m != Ainv.n || A.m != Ainv.m) return false;
        
        int n = A.m;
        Matrix LU = A; // Copy A
        Matrix p(n, 1);
        if (!lu_piv(LU, p)) return false;

        zeros(Ainv);
        for (int j = 0; j < n; j++) {
            // Solve A * x_j = e_j
            Matrix ej(n, 1, 0.0);
            ej(j, 0) = 1.0;
            if (!solve_lu_piv(LU, p, ej)) return false;
            
            // Copy result to j-th column of Ainv
            for (int i = 0; i < n; i++) {
                Ainv(i, j) = ej(i, 0);
            }
        }
        return true;
    }

    bool chol(Matrix& A) {
        if (A.m != A.n) return false;
        int n = A.m;
        for (int k = 0; k < n; k++) {
            double alpha_kk = A(k, k);
            if (alpha_kk <= 0.0) return false; // Matrix not SPD
            alpha_kk = Maths::sqrt(alpha_kk);
            A(k, k) = alpha_kk;

            // Scale column below diagonal
            for (int i = k + 1; i < n; i++) {
                A(i, k) /= alpha_kk;
            }

            // Symmetric rank-1 update of trailing matrix
            if (k < n - 1) {
                Matrix a21 = A.slice(k + 1, n, k, k + 1);
                Matrix A22 = A.slice(k + 1, n, k + 1, n);
                syr_l(-1.0, a21, A22);
            }
        }
        // Zero out the upper triangle to be clean
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < j; i++) {
                A(i, j) = 0.0;
            }
        }
        return true;
    }

    bool solve_chol(const Matrix& L, Matrix& b) {
        if (L.m != L.n || b.m * b.n != L.m) return false;
        int n = L.m;

        // Forward substitution: Lz = b
        for (int i = 0; i < n; i++) {
            double& b_i = (b.m == 1) ? b(0, i) : b(i, 0);
            for (int j = 0; j < i; j++) {
                double b_j = (b.m == 1) ? b(0, j) : b(j, 0);
                b_i -= L(i, j) * b_j;
            }
            if (L(i, i) == 0.0) return false;
            b_i /= L(i, i);
        }

        // Backward substitution: L^T x = z (z is now in b)
        for (int i = n - 1; i >= 0; i--) {
            double& b_i = (b.m == 1) ? b(0, i) : b(i, 0);
            for (int j = i + 1; j < n; j++) {
                double b_j = (b.m == 1) ? b(0, j) : b(j, 0);
                b_i -= L(j, i) * b_j; // Using L(j, i) for L^T(i, j)
            }
            // L^T(i, i) is L(i, i)
            b_i /= L(i, i);
        }
        return true;
    }

    bool rref(Matrix& A, double tol) {
        int m = A.m;
        int n = A.n;
        int pivot_row = 0;

        for (int j = 0; j < n && pivot_row < m; j++) {
            // Find pivot in column j, starting from pivot_row
            int max_row = pivot_row;
            double max_val = (A(max_row, j) < 0) ? -A(max_row, j) : A(max_row, j);

            for (int i = pivot_row + 1; i < m; i++) {
                double val = (A(i, j) < 0) ? -A(i, j) : A(i, j);
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }

            if (max_val < tol) {
                // All entries in current column from pivot_row down are effectively zero
                for (int i = pivot_row; i < m; i++) A(i, j) = 0.0;
                continue;
            }

            // Swap rows
            swap_rows(A, pivot_row, max_row);

            // Scale pivot row
            double pivot_val = A(pivot_row, j);
            for (int k = j; k < n; k++) {
                A(pivot_row, k) /= pivot_val;
            }
            A(pivot_row, j) = 1.0; // Ensure it's exactly 1.0

            // Eliminate entries in pivot column
            for (int i = 0; i < m; i++) {
                if (i != pivot_row) {
                    double factor = A(i, j);
                    for (int k = j; k < n; k++) {
                        A(i, k) -= factor * A(pivot_row, k);
                    }
                    A(i, j) = 0.0; // Ensure it's exactly 0.0
                }
            }
            pivot_row++;
        }
        return true;
    }

    int rank(const Matrix& A, double tol) {
        Matrix tmp = A;
        rref(tmp, tol);
        int r = 0;
        for (int i = 0; i < tmp.m; i++) {
            bool row_is_zero = true;
            for (int j = 0; j < tmp.n; j++) {
                double val = (tmp(i, j) < 0) ? -tmp(i, j) : tmp(i, j);
                if (val > tol) {
                    row_is_zero = false;
                    break;
                }
            }
            if (!row_is_zero) r++;
        }
        return r;
    }

    bool is_linearly_independent(const Matrix& A, double tol) {
        // iff rank(A) == A.n (number of columns)
        return rank(A, tol) == A.n;
    }

    bool is_spanning(const Matrix& A, double tol) {
        // iff rank(A) == A.m (dimension of target space R^m)
        return rank(A, tol) == A.m;
    }

    bool is_basis(const Matrix& A, double tol) {
        // iff A is square and rank(A) == A.n
        return (A.m == A.n) && (rank(A, tol) == A.n);
    }

    void print_matrix(const Matrix& A) {
        for (int i = 0; i < A.m; i++) {
            std::cout << "[ ";
            for (int j = 0; j < A.n; j++) {
                std::cout << std::fixed << std::setprecision(2) << std::setw(8) << A(i, j) << " ";
            }
            std::cout << " ]\n";
        }
    }
}