#include "Maths.hpp"

/*
 * Implemented by Shoko on 2026/03/06
 */
namespace Maths {

    double sqrt(double value) {
        if (value < 0.0) {
            return -1.0;
        }
        if (value == 0.0) {
            return 0.0;
        }

        double res = value;
        const double precision = 1e-15;

        for (int i = 0; i < 100; ++i) {
            double next = 0.5 * (res + value / res);
            double diff = (next > res) ? (next - res) : (res - next);
            if (diff < precision) {
                break;
            }
            res = next;
        }
        return res;
    }
} // namespace Maths
