#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Vector.hpp"

/**
 * Basic header of matrices
 * Currently in row major, but i might implement column-major later.
 * Created the 15/02/2026 by Shoko_ofi
 */

struct Matrix {
  int m;
  int n;

  Vector data;
  Matrix(int rows, int cols);

  double &operator()(int i, int j);
  const double &operator()(int i, int j) const;
};

#endif // MATRIX_HPP