#include "Laff.hpp"
#include "Matrix.hpp"
#include <iomanip>
#include <iostream>

/*
 * Test suite for LAFF routines
 */

void print_res(const char* label, bool ok, const Matrix& res, double expected = 0.0, bool is_scalar = false) {
    std::cout << (ok ? "[OK] " : "[FAIL] ") << label << "\n";
    if (ok) {
        if (is_scalar) {
            std::cout << "Value: " << res.data[0] << " (Expected: " << expected << ")\n";
        } else {
            std::cout << "Result: [ ";
            int size = res.m * res.n;
            for (int i = 0; i < size; ++i) {
                std::cout << res.data[i] << (i == size - 1 ? "" : ", ");
            }
            std::cout << " ]\n";
        }
    }
    std::cout << "----------------------------\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n=== LAFF Testing Suite ===\n";

    Matrix x(2, 1);
    x.data[0] = 4.0;
    x.data[1] = -3.0;
    Matrix y(2, 1);
    
    bool s1 = laff::copy(x, y);
    print_res("Vector Copy", s1, y);

    Matrix v_scal(2, 1);
    v_scal.data[0] = -1.0;
    v_scal.data[1] = 2.0;
    
    bool s2 = laff::scal(3.0, v_scal);
    print_res("Scaling (3.0 * [-1, 2])", s2, v_scal);

    Matrix x_axpy(2, 1);
    x_axpy.data[0] = -1.0;
    x_axpy.data[1] = 2.0;
    Matrix y_axpy(2, 1);
    y_axpy.data[0] = -3.0;
    y_axpy.data[1] = -2.0;

    bool s3 = laff::axpy(1.0, x_axpy, y_axpy);
    print_res("AXPY (Addition)", s3, y_axpy);

    Matrix v1(4, 1), v2(4, 1);
    v1.data[0] = 2.0; v1.data[1] = 5.0; v1.data[2] = -6.0; v1.data[3] = 1.0;
    for (int i = 0; i < 4; ++i) v2.data[i] = 1.0;

    double dot_out = 0.0;
    bool s4 = laff::dot(v1, v2, dot_out);
    Matrix tmp_dot(1, 1);
    tmp_dot.data[0] = dot_out;
    print_res("Dot Product", s4, tmp_dot, 2.0, true);

    Matrix v_norm(3, 1);
    v_norm.data[0] = 1.0; v_norm.data[1] = -2.0; v_norm.data[2] = 2.0;
    
    double norm_out = 0.0;
    bool s5 = laff::norm2(v_norm, norm_out);
    Matrix tmp_norm(1, 1);
    tmp_norm.data[0] = norm_out;
    print_res("L2 Norm", s5, tmp_norm, 3.0, true);

    Matrix small(2, 1), large(3, 1);
    bool error = laff::axpy(1.0, small, large);
    std::cout << (!error ? "[OK] " : "[FAIL] ") << "Size mismatch handling\n";
    std::cout << "----------------------------\n";

    // --- Matrix Special Tests ---
    std::cout << "\n=== Matrix Special Tests ===\n";

    Matrix A_zeros(2, 2, 5.0);
    laff::zeros(A_zeros);
    print_res("Matrix Zeros (2x2)", true, A_zeros);

    Matrix A_ident(3, 3);
    laff::identity(A_ident);
    print_res("Matrix Identity (3x3)", true, A_ident);

    Matrix v_diag(3, 1);
    v_diag.data[0] = 1.0; v_diag.data[1] = 2.0; v_diag.data[2] = 3.0;
    Matrix A_diag(3, 3);
    laff::diag(v_diag, A_diag);
    print_res("Matrix Diagonal (1,2,3)", true, A_diag);

    Matrix A_src(2, 3);
    A_src(0,0)=1.0; A_src(0,1)=2.0; A_src(0,2)=3.0;
    A_src(1,0)=4.0; A_src(1,1)=5.0; A_src(1,2)=6.0;
    Matrix B_trans(3, 2);
    laff::transpose(A_src, B_trans);
    print_res("Matrix Transpose (2x3 -> 3x2)", true, B_trans);

    Matrix A_tril(3, 3, 1.0);
    laff::tril(A_tril);
    print_res("Tril (Lower Triangular)", true, A_tril);

    Matrix A_triu(3, 3, 1.0);
    laff::triu(A_triu);
    print_res("Triu (Upper Triangular)", true, A_triu);

    std::cout << "Tests completed.\n\n";
    return 0;
}
