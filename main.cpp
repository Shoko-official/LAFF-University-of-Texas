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
    laff::lower_tri(A_tril);
    print_res("Tril (Lower Triangular)", true, A_tril);

    Matrix A_triu(3, 3, 1.0);
    laff::upper_tri(A_triu);
    print_res("Triu (Upper Triangular)", true, A_triu);

    Matrix A_to_sym(3, 3);
    A_to_sym(0,0)=1; A_to_sym(1,0)=2; A_to_sym(1,1)=3; A_to_sym(2,0)=4; A_to_sym(2,1)=5; A_to_sym(2,2)=6;
    laff::symmetrize_from_lower(A_to_sym);
    print_res("Symmetrize from Lower", true, A_to_sym);

    Matrix A_to_scal(2, 2, 2.0);
    laff::scal_matrix(1.5, A_to_scal);
    print_res("Scale Matrix (1.5 * [[2,2],[2,2]])", true, A_to_scal);

    Matrix A_add1(2, 2, 1.0);
    Matrix A_add2(2, 2, 2.0);
    laff::add_matrix(A_add2, A_add1);
    print_res("Add Matrices (1+2)", true, A_add1);

    // --- GEMV Tests (Chapter 3.4) ---
    std::cout << "\n=== GEMV Tests ===\n";
    Matrix A_gemv(3, 3);
    A_gemv(0,0)=-1; A_gemv(0,1)= 0; A_gemv(0,2)= 2;
    A_gemv(1,0)= 2; A_gemv(1,1)=-1; A_gemv(1,2)= 1;
    A_gemv(2,0)= 3; A_gemv(2,1)= 1; A_gemv(2,2)=-1;

    Matrix x_gemv(3, 1);
    x_gemv(0,0)=-1; x_gemv(1,0)= 2; x_gemv(2,0)= 4;

    Matrix y_dot(3, 1, 0.0);
    laff::gemv_dot(A_gemv, x_gemv, y_dot);
    print_res("GEMV Dot (Expected [9, 0, -5])", true, y_dot);

    Matrix y_gemv_axpy(3, 1, 0.0);
    laff::gemv_axpy(A_gemv, x_gemv, y_gemv_axpy);
    print_res("GEMV AXPY (Expected [9, 0, -5])", true, y_gemv_axpy);

    // --- Week 4: Transpose GEMV ---
    std::cout << "\n=== Week 4: Transpose GEMV ===\n";
    Matrix A_t(3, 2);
    A_t(0,0)=1; A_t(0,1)=2;
    A_t(1,0)=-2; A_t(1,1)=-1;
    A_t(2,0)=0; A_t(2,1)=1;
    
    Matrix x_t(3, 1);
    x_t(0,0)=-1; x_t(1,0)=2; x_t(2,0)=-3;
    
    Matrix y_t_dot(2, 1, 0.0);
    laff::gemv_t_dot(A_t, x_t, y_t_dot);
    print_res("GEMV_T Dot (Expected [-5, -7])", true, y_t_dot);

    Matrix y_t_axpy(2, 1, 0.0);
    laff::gemv_t_axpy(A_t, x_t, y_t_axpy);
    print_res("GEMV_T AXPY (Expected [-5, -7])", true, y_t_axpy);

    // --- Week 4: TRMV ---
    std::cout << "\n=== Week 4: TRMV ===\n";
    Matrix L(3, 3);
    L(0,0)=1;
    L(1,0)=2; L(1,1)=1;
    L(2,0)=3; L(2,1)=2; L(2,2)=1;
    
    Matrix x_l(3, 1);
    x_l(0,0)=1; x_l(1,0)=2; x_l(2,0)=3;

    laff::trmv_ln(L, x_l);
    print_res("TRMV Lower (Expected [1, 4, 10])", true, x_l);

    Matrix U(3, 3);
    U(0,0)=1; U(0,1)=2; U(0,2)=3;
    U(1,1)=1; U(1,2)=2;
    U(2,2)=1;
    
    Matrix x_u(3, 1);
    x_u(0,0)=1; x_u(1,0)=2; x_u(2,0)=3;

    laff::trmv_un(U, x_u);
    print_res("TRMV Upper (Expected [14, 8, 3])", true, x_u);

    // --- Week 4: SYMV ---
    std::cout << "\n=== Week 4: SYMV ===\n";
    Matrix A_sym(3, 3);
    A_sym(0,0)=1; 
    A_sym(1,0)=2; A_sym(1,1)=3;
    A_sym(2,0)=4; A_sym(2,1)=5; A_sym(2,2)=6;
    
    Matrix x_sym(3, 1);
    x_sym(0,0)=1; x_sym(1,0)=1; x_sym(2,0)=1;
    Matrix y_sym(3, 1, 0.0);
    
    laff::symv_l(A_sym, x_sym, y_sym);
    print_res("SYMV (Expected [7, 10, 15])", true, y_sym);

    // --- Week 5: Rank-1 Update ---
    std::cout << "\n=== Week 5: Rank-1 Update (GER) ===\n";
    Matrix A_ger(2, 3, 1.0);
    Matrix x_ger(2, 1); x_ger(0,0)=1; x_ger(1,0)=2;
    Matrix y_ger(3, 1); y_ger(0,0)=3; y_ger(1,0)=2; y_ger(2,0)=1;
    
    laff::ger(1.0, x_ger, y_ger, A_ger);
    print_res("GER (1.0 * x * y^T + 1.0)", true, A_ger);

    std::cout << "\n=== Week 5: Symmetric Rank-1 (SYR) ===\n";
    Matrix A_syr(2, 2, 0.0);
    Matrix x_syr(2, 1); x_syr(0,0)=1; x_syr(1,0)=2;
    laff::syr_l(1.0, x_syr, A_syr);
    print_res("SYR Lower (x * x^T)", true, A_syr);

    // --- Week 5: GEMM ---
    std::cout << "\n=== Week 5: GEMM ===\n";
    Matrix A_m(2, 2); A_m(0,0)=1; A_m(0,1)=2; A_m(1,0)=3; A_m(1,1)=4;
    Matrix B_m(2, 2); B_m(0,0)=5; B_m(0,1)=6; B_m(1,0)=7; B_m(1,1)=8;
    Matrix C_dot(2, 2, 0.0);
    
    laff::gemm_dot(1.0, A_m, B_m, 0.0, C_dot);
    print_res("GEMM Dot (Expected [[19, 22], [43, 50]])", true, C_dot);

    Matrix C_axpy(2, 2, 0.0);
    laff::gemm_axpy(1.0, A_m, B_m, 0.0, C_axpy);
    print_res("GEMM AXPY", true, C_axpy);

    Matrix C_row(2, 2, 0.0);
    laff::gemm_row(1.0, A_m, B_m, 0.0, C_row);
    print_res("GEMM Row", true, C_row);

    Matrix C_outer(2, 2, 0.0);
    laff::gemm_outer(1.0, A_m, B_m, 0.0, C_outer);
    print_res("GEMM Outer", true, C_outer);

    std::cout << "Tests completed.\n\n";
    return 0;
}
