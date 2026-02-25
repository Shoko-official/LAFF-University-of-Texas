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
}

#endif // LAFF_HPP
