#ifndef LAFF_HPP
#define LAFF_HPP

#include "Matrix.hpp"

/*
 * Linear Algebra: Foundations to Frontiers (LAFF) routines.
 * Regarding the UT Austin Course:
 * https://www.cs.utexas.edu/~flame/laff/laff/Notes_edX.pdf#page=31
 * 
 * Created by Shoko on 2026/03/06
 */
namespace laff {

    /**
     * Copies the contents of vector x into vector y of the same size
     * Handles all combinations of row and column vectors.
     * @param x Source vector
     * @param y Destination vector
     * @return True: copy succeeded, false if dimensions are invalid or incompatible
     */
    bool copy(const Matrix& x, Matrix& y);

    /**
     * Scales vector x by scalar alpha.
     * Maths:  x := ax
     * @param alpha Scaling factor
     * @param x Vector to be scaled in place
     * @return True: scaling succeeded, False: if x is not a valid vector
     */
    bool scal(double alpha, Matrix& x);

    /**
     * Computes y := alpha * x + y.
     * See page 38
     * @param alpha Scaling factor
     * @param x Source vector
     * @param y Destination vector (updated in place)
     * @return True: axpy operation succeeded, false if dimensions are invalid or incompatible
     */
    bool axpy(double alpha, const Matrix& x, Matrix& y);

    /**
     * Computes the dot product of two vectors: alpha := x^T y.
     * See page 40
     * @param x Source vector
     * @param y Destination vector
     * @param alpha Resulting scalar
     * @return True : dot operation succeeded, false if dimensions are invalid of incompatibl
     */
    bool dot(const Matrix& x, const Matrix& y, double& alpha);

     /**
     * Computes the Euclidean length (two-norm) of a vector: alpha := ||x||_2.
     * See page 43
     * @param x Source vector
     * @param alpha Resulting scalar
     * @return True if successful, false if x is not a vector.
     */
    bool norm2(const Matrix& x, double& alpha);

    /**
     * Fills the matrix A with zeros.
     * @param A Matrix to be zeroed
     * @return True if successful
     */
    bool zeros(Matrix& A);

    /**
     * Transforms matrix A into an identity matrix.
     * Matrix must be square (m == n).
     * @param A Matrix to be transformed
     * @return True if successful, false if A is not square
     */
    bool identity(Matrix& A);

    /**
     * Fills the diagonal of matrix A with elements from vector x.
     * A is set to zero elsewhere.
     * @param x Vector containing diagonal elements
     * @param A Matrix to be updated
     * @return True if successful, false if dimensions are incompatible
     */
    bool diag(const Matrix& x, Matrix& A);

    /**
     * Transposes matrix A and stores the result in B.
     * B must have dimensions A.n x A.m.
     * @param A Source matrix
     * @param B Destination matrix
     * @return True if successful, false if dimensions are incompatible
     */
    bool transpose(const Matrix& A, Matrix& B);

    /**
     * Sets all elements above the main diagonal of A to zero.
     * @param A Matrix to be processed in-place
     * @return True if successful
     */
    bool lower_tri(Matrix& A);

    /**
     * Sets all elements below the main diagonal of A to zero.
     * @param A Matrix to be processed in-place
     * @return True if successful
     */
    bool upper_tri(Matrix& A);

    /**
     * Symmetrizes matrix A by copying the lower triangle to the upper triangle.
     * @param A Matrix to be symmetrized (must be square)
     * @return True if successful
     */
    bool symmetrize_from_lower(Matrix& A);

    /**
     * Symmetrizes matrix A by copying the upper triangle to the lower triangle.
     * @param A Matrix to be symmetrized (must be square)
     * @return True if successful
     */
    bool symmetrize_from_upper(Matrix& A);

    /**
     * Scales matrix A by scalar alpha.
     * @param alpha Scaling factor
     * @param A Matrix to be scaled in place
     * @return True if successful
     */
    bool scal_matrix(double alpha, Matrix& A);

    /**
     * Adds matrix B to matrix A (A := A + B).
     * @param B Source matrix
     * @param A Destination matrix (updated in place)
     * @return True if successful, false if dimensions are incompatible
     */
    bool add_matrix(const Matrix& B, Matrix& A);

    /**
     * General Matrix-Vector Multiplication: y := Ax + y.
     * Implemented using dot products (row by row).
     * @param A Matrix m x n
     * @param x Vector n x 1
     * @param y Vector m x 1 (updated in place)
     * @return True if successful
     */
    bool gemv_dot(const Matrix& A, const Matrix& x, Matrix& y);

    /**
     * General Matrix-Vector Multiplication: y := Ax + y.
     * Implemented using AXPY operations (column by column).
     * @param A Matrix m x n
     * @param x Vector n x 1
     * @param y Vector m x 1 (updated in place)
     */
    bool gemv_axpy(const Matrix& A, const Matrix& x, Matrix& y);

    /**
     * General Matrix-Vector Multiplication with Transpose: y := A^T x + y.
     * Implemented using dot products.
     * @param A Matrix m x n
     * @param x Vector m x 1
     * @param y Vector n x 1 (updated in place)
     * @return True if successful
     */
    bool gemv_t_dot(const Matrix& A, const Matrix& x, Matrix& y);

    /**
     * General Matrix-Vector Multiplication with Transpose: y := A^T x + y.
     * Implemented using AXPY operations.
     * @param A Matrix m x n
     * @param x Vector m x 1
     * @param y Vector n x 1 (updated in place)
     * @return True if successful
     */
    bool gemv_t_axpy(const Matrix& A, const Matrix& x, Matrix& y);

    /**
     * Triangular Matrix-Vector Multiplication: x := L x (Lower, Non-transpose).
     * @param L Lower triangular matrix n x n
     * @param x Vector n x 1 (updated in place)
     * @return True if successful
     */
    bool trmv_ln(const Matrix& L, Matrix& x);

    /**
     * Triangular Matrix-Vector Multiplication: x := U x (Upper, Non-transpose).
     * @param U Upper triangular matrix n x n
     * @param x Vector n x 1 (updated in place)
     * @return True if successful
     */
    bool trmv_un(const Matrix& U, Matrix& x);

    /**
     * Symmetric Matrix-Vector Multiplication: y := A x + y.
     * Where A is symmetric and only the lower triangle is used.
     * @param A Symmetric matrix n x n
     * @param x Vector n x 1
     * @param y Vector n x 1 (updated in place)
     * @return True if successful
     */
    bool symv_l(const Matrix& A, const Matrix& x, Matrix& y);

    /**
     * Rank-1 Update (Outer Product): A := alpha * x * y^T + A.
     * Maths: A := alpha * x * y^T + A
     * See Week 5, Unit 5.3.4 (Page 183/190)
     * @param alpha Scalar factor
     * @param x Vector m x 1
     * @param y Vector n x 1
     * @param A Matrix m x n (updated in place)
     * @return True if successful
     */
    bool ger(double alpha, const Matrix& x, const Matrix& y, Matrix& A);

    /**
     * Symmetric Rank-1 Update: A := alpha * x * x^T + A.
     * Updates only the lower triangle of A.
     * See Week 5, Unit 5.3.4 (Page 183)
     * @param alpha Scalar factor
     * @param x Vector n x 1
     * @param A Symmetric matrix n x n (updated in place)
     * @return True if successful
     */
    bool syr_l(double alpha, const Matrix& x, Matrix& A);

    /**
     * General Matrix-Matrix Multiplication: C := alpha * A * B + beta * C.
     * Default implementation (calls gemm_axpy).
     * See Week 5, Chapter 5.3 (Page 183)
     * @param alpha Scalar factor
     * @param A Matrix m x k
     * @param B Matrix k x n
     * @param beta Scalar factor (Beta * C)
     * @param C Matrix m x n (updated in place)
     * @return True if successful
     */
    bool gemm(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C);

    /**
     * Matrix-Matrix Multiplication variant using dot products.
     * Computes C(i,j) := alpha * row_i(A) * col_j(B) + beta * C(i,j)
     * See Week 5, Unit 5.3.1 (Page 184) - "Lots of Loops"
     */
    bool gemm_dot(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C);

    /**
     * Matrix-Matrix Multiplication variant using AXPY (Matrix-Vector).
     * Computes C via columns: col_j(C) := alpha * A * col_j(B) + beta * col_j(C)
     * See Week 5, Unit 5.3.2 (Page 186/187) - "Matrix-Matrix Multiplication by Columns"
     * Corresponds to GEMM_UNB_VAR1 in FLAME notation.
     */
    bool gemm_axpy(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C);

    /**
     * Matrix-Matrix Multiplication variant using Rows.
     * Computes C via rows: row_i(C) := alpha * row_i(A) * B + beta * row_i(C)
     * See Week 5, Unit 5.3.3 (Page 188/189) - "Matrix-Matrix Multiplication by Rows"
     * Corresponds to GEMM_UNB_VAR2 in FLAME notation.
     */
    bool gemm_row(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C);

    /**
     * Matrix-Matrix Multiplication variant using Outer Products (Rank-1 updates).
     * Computes C via sum of outer products: C := alpha * sum( col_k(A) * row_k(B) ) + beta * C
     * See Week 5, Unit 5.3.4 (Page 190/191) - "Rank-1 Updates"
     * Corresponds to GEMM_UNB_VAR3 in FLAME notation.
     */
    bool gemm_outer(double alpha, const Matrix& A, const Matrix& B, double beta, Matrix& C);

    /**
     * Triangular Solve: x := L^-1 * x (Lower Triangular, Non-transpose).
     * Solves Lx = b where L is lower triangular, stored in-place in A.
     * Assumes L is unit lower triangular if used with LU.
     * See Week 6, Unit 6.3.2 (Page 212)
     * @param L Lower triangular matrix n x n
     * @param x Vector n x 1 (updated in place)
     * @return True if successful
     */
    bool trsv_ln(const Matrix& L, Matrix& x);

    /**
     * Triangular Solve: x := U^-1 * x (Upper Triangular, Non-transpose).
     * Solves Ux = b where U is upper triangular, stored in-place in A.
     * See Week 6, Unit 6.3.3 (Page 214)
     * @param U Upper triangular matrix n x n
     * @param x Vector n x 1 (updated in place)
     * @return True if successful
     */
    bool trsv_un(const Matrix& U, Matrix& x);

    /**
     * LU Factorization (Unblocked): A := LU.
     * Computes the LU factorization of a square matrix A in-place.
     * L is unit lower triangular (diagonal of 1s not stored).
     * See Week 6, Unit 6.3.1 (Page 209)
     * @param A Matrix n x n (updated in place)
     * @return True if successful
     */
    bool lu_unb(Matrix& A);

    /**
     * Solves Ax = b using LU Factorization.
     * Computes LU factorization of A, then solves via forward and backward substitution.
     * A and b are updated in-place (A contains LU, b contains x).
     * @param A Matrix n x n
     * @param b Vector n x 1
     * @return True if successful
     */
    bool solve_lu(Matrix& A, Matrix& b);

    /**
     * Finds the index of the element with the largest absolute value in a vector.
     * See Week 7, Unit 7.2.4
     * @param x Source vector
     * @param index Resulting index (0-indexed)
     * @return True if successful
     */
    bool iamax(const Matrix& x, int& index);

    /**
     * Swaps two rows of a matrix.
     * @param A Matrix to be updated in-place
     * @param i First row index
     * @param j Second row index
     * @return True if successful
     */
    bool swap_rows(Matrix& A, int i, int j);

    /**
     * LU Factorization with Partial Pivoting: PA := LU.
     * Computes the LU factorization of a square matrix A with row swapping.
     * p stores the pivot indices (permutation information).
     * See Week 7, Unit 7.2.4 (Page 267)
     * @param A Matrix n x n (updated in place with L and U)
     * @param p Permutation vector n x 1
     * @return True if successful
     */
    bool lu_piv(Matrix& A, Matrix& p);

    /**
     * Applies the permutation recorded in p to vector b.
     * See Week 7, Unit 7.2.4 (Page 268)
     * @param p Permutation vector n x 1
     * @param b Vector n x 1 (updated in-place)
     * @return True if successful
     */
    bool apply_piv(const Matrix& p, Matrix& b);

    /**
     * Solves Ax = b using LU Factorization with Partial Pivoting.
     * Updates A (LU) and b (x) in-place.
     * @param A Matrix n x n
     * @param p Permutation vector n x 1
     * @param b Vector n x 1
     * @return True if successful
     */
    bool solve_lu_piv(Matrix& A, const Matrix& p, Matrix& b);

    /**
     * Computes the inverse of a square matrix A using LU with pivoting.
     * @param A Source matrix n x n
     * @param Ainv Resulting inverse matrix n x n
     * @return True if successful (matrix is non-singular)
     */
    bool inv(const Matrix& A, Matrix& Ainv);

    /**
     * Cholesky Factorization: A := LL^T.
     * Computes the Cholesky factorization of a symmetric positive definite matrix A.
     * Overwrites the lower triangular part with L.
     * See Week 8, Unit 8.4.2 (Page 320)
     * @param A Matrix n x n (updated in place)
     * @return True if successful
     */
    bool chol(Matrix& A);

    /**
     * Solves Ax = b for a symmetric positive definite matrix A using its Cholesky factor L.
     * Assumes A has been factored using chol().
     * @param L Factored matrix n x n (only lower triangle used)
     * @param b Vector n x 1 (updated in place with x)
     * @return True if successful
     */
    bool solve_chol(const Matrix& L, Matrix& b);

    /**
     * @brief Computes the Reduced Row Echelon Form (RREF) of matrix A.
     * Uses Gauss-Jordan elimination with partial pivoting.
     * See Week 9/10 concepts.
     * @param A Matrix m x n (updated in-place)
     * @param tol Tolerance for zero detection (default 1e-10)
     * @return True if successful
     */
    bool rref(Matrix& A, double tol = 1e-10);

    /**
     * @brief Computes the rank of matrix A.
     * Rank is the number of non-zero rows in RREF.
     * @param A Matrix m x n
     * @param tol Tolerance for zero detection
     * @return The rank of the matrix
     */
    int rank(const Matrix& A, double tol = 1e-10);

    /**
     * @brief Checks if the columns of matrix A are linearly independent.
     * iff rank(A) == A.n
     */
    bool is_linearly_independent(const Matrix& A, double tol = 1e-10);

    /**
     * @brief Checks if the columns of A span the space R^m.
     * iff rank(A) == A.m
     */
    bool is_spanning(const Matrix& A, double tol = 1e-10);

    /**
     * @brief Checks if the columns of A form a basis for R^m.
     * iff A is square and rank(A) == A.n
     */
    bool is_basis(const Matrix& A, double tol = 1e-10);

    /**
     * @brief Prints the matrix A to standard output.
     * @param A Matrix to print
     */
    void print_matrix(const Matrix& A);
}

#endif // LAFF_HPP
