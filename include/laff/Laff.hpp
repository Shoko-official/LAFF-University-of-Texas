#ifndef LAFF_HPP
#define LAFF_HPP

#include "laff/Matrix.hpp"
#include "laff/Vector.hpp"
#include "laff/blas/VectorOperations.hpp"
#include "laff/blas/MatrixVectorOps.hpp"
#include "laff/blas/MatrixMatrixMult.hpp"
#include "laff/linear_systems/LU.hpp"
#include "laff/cholesky/Cholesky.hpp"
#include "laff/vector_spaces/VectorSpaces.hpp"
#include "laff/utility/Maths.hpp"
#include "laff/utility/Printing.hpp"

namespace laff {
    // Core utility routines
    bool fill(Matrix& A, double val);

    bool zeros(Matrix& A);
    bool identity(Matrix& A);


    bool diag(const Matrix& x, Matrix& A);
    bool transpose(const Matrix& A, Matrix& B);
    bool lower_tri(Matrix& A);
    bool upper_tri(Matrix& A);
    bool symmetrize_from_lower(Matrix& A);
    bool symmetrize_from_upper(Matrix& A);
    bool scal_matrix(double alpha, Matrix& A);
    bool add_matrix(const Matrix& B, Matrix& A);
}

#endif // LAFF_HPP
