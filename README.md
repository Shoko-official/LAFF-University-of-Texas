# LAFF-C++ | Linear Algebra Engine from First Principles

This repository contains my personal implementation of a C++ Linear Algebra library, developed while following the **LAFF (Linear Algebra: Foundations to Frontiers)** course [^1].

While I am using these routines as the backbone for a "LLM from scratch" project, this library stands on its own as a rigorous exploration of BLAS-like operations and matrix factorizations. The goal was simple but ambitious: stop treating `numpy` or `torch.matmul` as black boxes and understand the memory patterns and arithmetic required to move data efficiently.

## Why build this?

Most ML students start at the "Layer" level. I wanted to start at the "Memory" level.<br> 
Implementing $C = AB$ or LU decomposition in C++ forces to confront issues that high-level Python libraries hide:
* **Row-major vs Column-major** storage and its impact on cache hits.
* **Numerical stability** during pivoting for linear systems.
* **Pointer arithmetic** and avoiding memory leaks in dynamic matrix allocation.

## Current Feature Set

### 1. Fundamental BLAS Operations
I've implemented the standard hierarchy of operations:
* **Level 1**: Vector scaling, addition (AXPY), and dot products.
* **Level 2**: Matrix-vector multiplication (GEMV) using different algorithmic variants (dot-product vs. axis-based).
* **Level 3**: Matrix-Matrix multiplication (GEMM). This is the core engine for any future Neural Network implementation.

### 2. Linear Systems & Factorizations
To go beyond simple arithmetic, I've integrated solvers for more complex problems:
* **LU Decomposition**: Implemented with partial pivoting to handle systems of equations without numerical collapse.
* **Cholesky Factorization**: Specifically for symmetric positive definite matrices, useful for certain optimization algorithms.
* **Vector Spaces**: Functions to compute rank, RREF, and basis properties.

## Project Structure

* [/include](include/laff/): Header files defining the API. I opted for a clean separation between the [Matrix class](include/laff/Matrix.hpp) and [Vector class](include/laff/Vector.hpp) to keep the library modular and intuitive.
* [/src](src/laff/): Implementation of the computational kernels. You will find the core logic for [Matrix-Matrix Multiplication](src/laff/blas/MatrixMatrixMult.cpp) and the [LU Factorization](src/laff/linear_systems/LU.cpp) algorithms here.
* [/tests](tests/): A [comprehensive suite](tests/main.cpp) of manual tests to verify the correctness of every routine (from Level 1 BLAS to Cholesky) against known mathematical results.

## Technical Takeaways

One of the main challenges was managing the life cycle of large matrices. I focused on building a robust `Matrix` class that handles its own memory, ensuring that even when we chain operations (like in a forward pass of a model), we don't end up with dangling pointers or overhead.

## How to Compile
> [!TIP]
>I kept it lightweight. No heavy dependencies required.
>```bash
># Compile everything with the test runner
>g++ -Iinclude src/laff/**/*.cpp tests/main.cpp -o linalg_test
>./linalg_test
>```

## Contributions

Contributions

This is a solo research project, but I welcome technical discussions, especially regarding:
- Optimization: Ideas on Tiling or SIMD instructions to further speed up the gemm kernels.
- Numerical Analysis: Edge cases in partial pivoting or Cholesky stability.

If you have any questions, suggestions for optimization, or if you are a recruiter/professor from **Sorbonne Université**, please feel free to [open a technical issue here](https://github.com/shoko-official/laff-university-of-texas/issues/new).

---
[^1]: [certification](https://courses.edx.org/certificates/b4b26513557c47d2b199ae23c6905820)

*Authored by Shoko_official. I recommand you have a look at this other repository : [LLM From Scratch](https://github.com/Shoko-official/LLM-From-Abs-Scratch)*
