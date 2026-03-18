#include "laff/Laff.hpp"

namespace laff {
    bool zeros(Matrix& A) {
        for (int i = 0; i < A.m * A.n; i++) A.data[i] = 0.0;
        return true;
    }
    // ... (only core routines kept here)
}