#include "Laff.hpp"
#include "Utility/Maths.hpp"

/*
* Implemented by Shoko on 2026-03-06
*/

namespace laff {
    bool copy(const Matrix& x, Matrix& y) {
        if ((x.m != 1 && x.n != 1) || (y.m != 1 && y.n != 1)) {
            return false;
        }

        if (x.m * x.n != y.m * y.n) {
            return false;
        }

        int size = x.m * x.n;
        for (int i = 0; i < size; i++) {
            y.data[i] = x.data[i];
        }
           return true;
    }


    bool scal(double alpha, Matrix& x){
        if (x.m != 1 && x.n != 1) {
            return false;
        }

        int size = x.m * x.n;
        /* Choosed to bypass 2D operatior, to iterate directly over 1D 
        contigous memory array
        It avoid arithmetic overhead of indexes calculated by iteration, and guarantees sequential memory access.
        */
        for (int i = 0; i < size; i++) {
            x.data[i] *= alpha;
        }

        return true;
    }

    bool axpy(double alpha, const Matrix& x, Matrix& y) {
        if ((x.m != 1 && x.n != 1) || (y.m != 1 && y.n != 1)) {
            return false;
        }
        if (x.m * x.n != y.m * y.n) {
            return false;
        }
        int size = x.m * x.n;

        /* Choosed to bypass 2D operatior, to iterate directly over 1D 
        contigous memory array
        It avoid arithmetic overhead of indexes calculated by iteration, and guarantees sequential memory access.
        */
        for (int i = 0; i < size; i++) {
            // y := alpha * x + y
            y.data[i] = alpha * x.data[i] + y.data[i];
        }
        return true;
    }

    bool dot(const Matrix& x, const Matrix& y, double& alpha){
        if ((x.m != 1 && x.n != 1) || (y.m != 1 && y.n != 1)) {
            return false;
        }
        if (x.m * x.n != y.m * y.n) {
            return false;
        }
        int size = x.m * x.n;
        alpha = 0.0;
        
        for (int i = 0; i < size; i++){
            alpha += x.data[i] * y.data[i];
        }
        return true;
    }

    bool norm2(const Matrix& x, double& alpha) {
        if (x.m != 1 && x.n != 1) {
            return false;
        }

        if (!dot(x, x, alpha)) {
            return false;
        }

        alpha = Maths::sqrt(alpha);
        return true;
    }

    bool zeros(Matrix& A) {
        int size = A.m * A.n;
        /* Choosed to bypass 2D operatior, to iterate directly over 1D 
        contigous memory array as it avoids arithmetic overhead of indexes calculated by iteration
		and guarantees sequential memory access.
        */
        for (int i = 0; i < size; i++) {
            A.data[i] = 0.0;
        }
        return true;
    }

    bool identity(Matrix& A) {
        if (A.m != A.n) {
            return false;
        }

        zeros(A);

        /* Set the diagonal elements to 1.0 */
        for (int i = 0; i < A.m; i++) {
            A(i, i) = 1.0;
        }
        return true;
    }
}
