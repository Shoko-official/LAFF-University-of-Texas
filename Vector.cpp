#include "Vector.hpp"

/**
 * Basic implementation of vectors
 * Created the 15/02/2026 by Shoko_ofi
 */

Vector::Vector(int size, double val) {
    length = size;
    capacity = (size > 0) ? size : 4;
    data = new double[capacity];
    for (int i = 0; i < length; i++) {
        data[i] = val;
    }
}

Vector::~Vector() {
    delete[] data;
}

Vector::Vector(const Vector& other) : length(other.length), capacity(other.capacity) {
    data = new double[capacity];
    for (int i = 0; i < length; ++i) {
        data[i] = other.data[i];
    }
}

Vector& Vector::operator=(const Vector& other) {
    if (this != &other) {
        delete[] data;
        length = other.length;
        capacity = other.capacity;
        data = new double[capacity];
        for (int i = 0; i < length; ++i) {
            data[i] = other.data[i];
        }
    } return *this;
}
