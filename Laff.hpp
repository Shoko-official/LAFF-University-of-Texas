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
}

#endif // LAFF_HPP