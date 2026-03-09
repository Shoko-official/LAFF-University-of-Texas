#ifndef LU_HPP
#define LU_HPP

#include "../Matrix.hpp"

namespace laff {
    bool lu_unb(Matrix& A);
    bool solve_lu(Matrix& A, Matrix& b);
    bool iamax(const Matrix& x, int& index);
    bool swap_rows(Matrix& A, int i, int j);
    bool lu_piv(Matrix& A, Matrix& p);
    bool apply_piv(const Matrix& p, Matrix& b);
    bool solve_lu_piv(Matrix& A, const Matrix& p, Matrix& b);
    bool inv(const Matrix& A, Matrix& Ainv);
}

#endif // LU_HPP
