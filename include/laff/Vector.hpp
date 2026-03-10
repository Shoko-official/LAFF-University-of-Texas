#ifndef VECTOR_HPP
#define VECTOR_HPP

/**
 * Basic header of vectors
 * Created the 2026/02/15 by Shoko_ofi
 * Editted the 2026/03/06 by Shoko_ofi
 */

struct Vector {
    double* data;
    int length;
    int capacity;
    int count;
    bool owns_memory;

    Vector(int size = 0, double val = 0.0);
    Vector(double* prt, int size, int stride);
    ~Vector();

    Vector(const Vector& other);
    Vector& operator=(const Vector& other);

    double& operator[](int index);

    const double& operator[](int index) const;

    Vector slice(int start, int end);
};

#endif // VECTOR_HPP