#ifndef LAFF_HPP
#define LAFF_HPP

#include "Matrix.hpp"
#include "Vector.hpp"
#include "blas/VectorOperations.hpp"
#include "blas/MatrixVectorOps.hpp"
#include "blas/MatrixMatrixMult.hpp"
#include "linear_systems/LU.hpp"
#include "cholesky/Cholesky.hpp"
#include "vector_spaces/VectorSpaces.hpp"
#include "utility/Maths.hpp"
#include "utility/Printing.hpp"

namespace laff {
    // Core utility routines
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
