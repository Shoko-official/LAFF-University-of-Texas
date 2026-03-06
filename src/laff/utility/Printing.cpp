#include "laff/Laff.hpp"
#include <iostream>
#include <iomanip>

namespace laff {
    void print_matrix(const Matrix& A) {
        for (int i = 0; i < A.m; i++) {
            std::cout << "[ ";
            for (int j = 0; j < A.n; j++) {
                std::cout << std::setw(10) << A(i, j) << " ";
            }
            std::cout << "]" << std::endl;
        }
    }
}
