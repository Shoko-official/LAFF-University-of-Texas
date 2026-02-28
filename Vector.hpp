#ifndef VECTOR_HPP
#define VECTOR_HPP

/**
 * Basic header of vectors
 * Created the 15/02/2026 by Shoko_ofi
 */

struct Vector {
    double* data;
    int length;
    int capacity;

    Vector(int size = 0, double val = 0.0);
    ~Vector();

    Vector(const Vector& other);
    Vector& operator=(const Vector& other);

    double& operator[](int index) {
        return data[index];
    }

    const double& operator[](int index) const {
        return data[index];
    }
};

#endif // VECTOR_HPP