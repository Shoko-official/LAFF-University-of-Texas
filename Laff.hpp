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
}

#endif // LAFF_HPP