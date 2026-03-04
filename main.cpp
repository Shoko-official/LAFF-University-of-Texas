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

    // --- Week 6: LU and Solvers ---
    std::cout << "\n=== Week 6: LU and Solvers ===\n";
    
    // Test TRSV Forward
    Matrix L_sol(3, 3, 0.0);
    L_sol(0,0)=1; 
    L_sol(1,0)=2; L_sol(1,1)=1;
    L_sol(2,0)=-1; L_sol(2,1)=2; L_sol(2,2)=1;
    Matrix b_l(3, 1); b_l(0,0)=1; b_l(1,0)=4; b_l(2,0)=9;
    laff::trsv_ln(L_sol, b_l);
    print_res("TRSV Forward (Expected [1, 2, 6])", true, b_l);

    // Test TRSV Backward
    Matrix U_sol(3, 3, 0.0);
    U_sol(0,0)=2; U_sol(0,1)=4; U_sol(0,2)=-2;
    U_sol(1,1)=-10; U_sol(1,2)=10;
    U_sol(2,2)=-8;
    Matrix b_u(3, 1); b_u(0,0)=-10; b_u(1,0)=40; b_u(2,0)=-16;
    laff::trsv_un(U_sol, b_u);
    print_res("TRSV Backward (Expected [1, -2, 2])", true, b_u);

    // Test LU and Full Solver
    Matrix A_solve(3, 3);
    A_solve(0,0)=2; A_solve(0,1)=4; A_solve(0,2)=-2;
    A_solve(1,0)=4; A_solve(1,1)=-2; A_solve(1,2)=6;
    A_solve(2,0)=6; A_solve(2,1)=-4; A_solve(2,2)=2;
    Matrix b_solve(3, 1); b_solve(0,0)=-10; b_solve(1,0)=20; b_solve(2,0)=18;
    
    laff::solve_lu(A_solve, b_solve);
    print_res("Full Solve_LU (Expected [1, -2, 2])", true, b_solve);

    // --- Week 7: Partial Pivoting and Inversion ---
    std::cout << "\n=== Week 7: Partial Pivoting and Inversion ===\n";

    // Matrix requiring pivot (A(0,0)=0)
    Matrix A_swap(2, 2);
    A_swap(0,0)=0; A_swap(0,1)=1;
    A_swap(1,0)=1; A_swap(1,1)=1;
    Matrix b_swap(2, 1);
    b_swap(0,0)=1; b_swap(1,0)=2; // System: 0x + 1y = 1, 1x + 1y = 2 -> Solution: x=1, y=1

    Matrix p_idx(2, 1);
    bool s_lup = laff::lu_piv(A_swap, p_idx);
    std::cout << (s_lup ? "[OK] " : "[FAIL] ") << "LU_PIV on matrix with zero pivot\n";

    laff::solve_lu_piv(A_swap, p_idx, b_swap);
    print_res("Solve_LU_PIV (Expected [1, 1])", true, b_swap);

    // Matrix Inversion
    Matrix A_inv_src(3, 3);
    A_inv_src(0,0)=2; A_inv_src(0,1)=1; A_inv_src(0,2)=0;
    A_inv_src(1,0)=1; A_inv_src(1,1)=2; A_inv_src(1,2)=1;
    A_inv_src(2,0)=0; A_inv_src(2,1)=1; A_inv_src(2,2)=2;
    
    Matrix A_inv_res(3, 3);
    bool s_inv = laff::inv(A_inv_src, A_inv_res);
    std::cout << (s_inv ? "[OK] " : "[FAIL] ") << "Matrix Inversion\n";
    
    // Check A * A_inv
    Matrix I_check(3, 3, 0.0);
    laff::gemm(1.0, A_inv_src, A_inv_res, 0.0, I_check);
    print_res("Identity Check (A * A^-1)", true, I_check);

    // --- Week 8: Cholesky Factorization ---
    std::cout << "\n=== Week 8: Cholesky Factorization ===\n";
    Matrix A_chol(2, 2);
    A_chol(0,0)=4; A_chol(0,1)=2;
    A_chol(1,0)=2; A_chol(1,1)=2;
    Matrix b_chol(2, 1);
    b_chol(0,0)=6; b_chol(1,0)=4;

    Matrix A_chol_copy = A_chol;
    bool s_chol = laff::chol(A_chol);
    std::cout << (s_chol ? "[OK] " : "[FAIL] ") << "Cholesky Factorization\n";
    print_res("L factor (Expected [[2,0], [1,1]])", s_chol, A_chol);

    laff::solve_chol(A_chol, b_chol);
    print_res("Solve_CHOL (Expected [1, 1])", s_chol, b_chol);

    std::cout << "\n=== Week 9: Vector Spaces ===\n";
    Matrix A_rref(3, 3);
    A_rref(0,0)=1; A_rref(0,1)=2; A_rref(0,2)=3;
    A_rref(1,0)=4; A_rref(1,1)=5; A_rref(1,2)=6;
    A_rref(2,0)=7; A_rref(2,1)=8; A_rref(2,2)=9;
    
    std::cout << "Original Matrix (3x3, Rank 2):\n";
    laff::print_matrix(A_rref);
    
    int r = laff::rank(A_rref);
    std::cout << "Rank (Expected 2): " << r << (r == 2 ? " [OK]" : " [FAIL]") << "\n";
    
    laff::rref(A_rref);
    std::cout << "RREF (Expected row 3 as zeros):\n";
    laff::print_matrix(A_rref);
    
    Matrix B_indep(3, 2);
    B_indep(0,0)=1; B_indep(0,1)=0;
    B_indep(1,0)=0; B_indep(1,1)=1;
    B_indep(2,0)=0; B_indep(2,1)=0;
    bool is_indep = laff::is_linearly_independent(B_indep);
    std::cout << "Is linearly independent (Expected 1): " << is_indep << (is_indep ? " [OK]" : " [FAIL]") << "\n";
    
    Matrix C_basis(2, 2);
    C_basis(0,0)=1; C_basis(0,1)=0;
    C_basis(1,0)=0; C_basis(1,1)=1;
    bool is_bas = laff::is_basis(C_basis);
    std::cout << "Is basis (Expected 1): " << is_bas << (is_bas ? " [OK]" : " [FAIL]") << "\n";

    std::cout << "Tests completed.\n\n";
    return 0;
}
