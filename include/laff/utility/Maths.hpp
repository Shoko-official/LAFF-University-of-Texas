#ifndef MATHS_HPP
#define MATHS_HPP

/*
 * Mathematical utility functions.
 * Created by Shoko on 2026/03/06
 */
namespace Maths {
/**
 * Computes the square root of a non-negative double using Newton's method.
 * See https://math.mit.edu/~stevenj/18.335/newton-sqrt.pdf
 * @param value The number to compute the square root of.
 * @return The square root, or -1.0 if the input is negative (edge case).
 */
double sqrt(double value);
} // namespace maths

#endif // MATHS_HPP