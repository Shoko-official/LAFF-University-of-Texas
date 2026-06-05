#ifndef VECTOR_HPP
#define VECTOR_HPP

/**
 * Basic header of vectors
 */

namespace laff {

template <typename T = double>
struct Vector {
    T* data;
    int length;
    int capacity;
    int count;
    bool owns_memory;

    Vector(int size = 0, T val = T(0)) : length(size), capacity(size > 0 ? size : 4), count(1), owns_memory(true) {
        data = new T[capacity];
        for (int i = 0; i < length; i++) {
            data[i] = val;
        }
    }

    Vector(T* ptr, int size, int stride) : data(ptr), length(size), capacity(size), count(stride), owns_memory(false) {}

    ~Vector() {
        if (owns_memory) {
            delete[] data;
        }
    }

    Vector(const Vector& other) : length(other.length), capacity(other.length), count(1), owns_memory(true) {
        data = new T[capacity];
        for (int i = 0; i < length; ++i) {
            data[i] = other.data[i * other.count];
        }
    }

    Vector& operator=(const Vector& other) {
        if (this != &other) {
            if (owns_memory) {
                delete[] data;
            }
            length = other.length;
            capacity = other.capacity;
            count = 1;
            owns_memory = true;
            data = new T[capacity];
            for (int i = 0; i < length; ++i) {
                data[i] = other.data[i * other.count];
            }
        }
        return *this;
    }

    T& operator[](int index) {
        return data[index * count];
    }

    const T& operator[](int index) const {
        return data[index * count];
    }

    Vector slice(int start, int end) {
        return Vector(data + (start * count), end - start, count);
    }
};

} // namespace laff

// For compatibility with global namespace:
using Vector = laff::Vector<double>;

#endif // VECTOR_HPP