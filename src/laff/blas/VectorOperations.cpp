#include "laff/Laff.hpp"
#include "laff/utility/Maths.hpp"

namespace laff {
    bool copy(const Matrix& x, Matrix& y) {
        if ((x.m != 1 && x.n != 1) || (y.m != 1 && y.n != 1)) return false;
        if (x.m * x.n != y.m * y.n) return false;
        for (int i = 0; i < x.m * x.n; i++) y.data[i] = x.data[i];
        return true;
    }

    bool scal(double alpha, Matrix& x) {
        if (x.m != 1 && x.n != 1) return false;
        for (int i = 0; i < x.m * x.n; i++) x.data[i] *= alpha;
        return true;
    }

    bool axpy(double alpha, const Matrix& x, Matrix& y) {
        if (x.m * x.n != y.m * y.n) return false;
        for (int i = 0; i < x.m * x.n; i++) y.data[i] += alpha * x.data[i];
        return true;
    }

    bool dot(const Matrix& x, const Matrix& y, double& alpha) {
        if (x.m * x.n != y.m * y.n) return false;
        alpha = 0.0;
        for (int i = 0; i < x.m * x.n; i++) alpha += x.data[i] * y.data[i];
        return true;
    }

    bool norm2(const Matrix& x, double& alpha) {
        double d;
        if (!dot(x, x, d)) return false;
        alpha = Maths::sqrt(d);
        return true;
    }
}
