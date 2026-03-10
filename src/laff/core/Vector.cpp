#include "laff/Vector.hpp"

/**
 * Basic implementation of vectors
 * Created the 2026/02/15 by Shoko_ofi
 * Editted the 2026/03/06 by Shoko_ofi
 */

Vector::Vector(int size, double val) : length(size), capacity(size > 0 ? size : 4), count(1), owns_memory(true) {
    data = new double[capacity];
    for (int i = 0; i < length; i++) {
        data[i] = val;
    }
}

Vector::Vector(double* ptr, int size, int stride) : data(ptr), length(size), capacity(size), count(stride), owns_memory(false) {}

Vector::~Vector() {
    if (owns_memory) {
        delete[] data;
    }
}

Vector::Vector(const Vector& other) : length(other.length), capacity(other.length), count(1), owns_memory(true) {
    data = new double[capacity];
    for (int i = 0; i < length; ++i) {
        data[i] = other.data[i * other.count];
    }
}

Vector& Vector::operator=(const Vector& other) {
    if (this != &other) {
        if (owns_memory) {
            delete[] data;
        }
        length = other.length;
        capacity = other.capacity;
        count = 1;
        owns_memory = true;
        data = new double[capacity];
        for (int i = 0; i < length; ++i) {
            data[i] = other.data[i * other.count];
        }
    } return *this;
}

double& Vector::operator[](int idx) {
    return data[idx * count];
}

const double& Vector::operator[](int idx) const {
    return data[idx * count];
}

Vector Vector::slice(int start, int end) {
    return Vector(data + (start * count), end -start, count);
}