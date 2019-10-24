#include <iostream>
#include <tuple>
#include "Base.h"
#include "Vector.h"
#include "Matrix.h"
#include <vector>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace mat_vec;

const double eps = 1e-5;

//// 1
TEST_CASE("Constructor vector") {
    Vector A(5, 0);
    REQUIRE(A.size() == 5);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i]) < eps);
    }
}

//// 2
TEST_CASE("Constructor copy vector") {
    Vector A(5, 3);
    Vector B(A);
    REQUIRE(B.size() == 5);
    for (int i = 0; i < B.size(); ++i) {
        REQUIRE(abs(B[i] - 3) < eps);
    }
}

//// 3
TEST_CASE("Operator = vector") {
    Vector A(5, 0);
    Vector B(1);
    B = A;
    REQUIRE(B.size() == 5);
    for (int i = 0; i < B.size(); ++i) {
        REQUIRE(abs(B[i]) < eps);
    }
}

//// 4
TEST_CASE("Size vector") {
    Vector A(5, 0);
    REQUIRE(A.size() == 5);
}

//// 5
TEST_CASE("Operator [] vector") {
    Vector A(5, 10);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 10) < eps);
    }
    double k = A[1];
    REQUIRE(abs(k - 10) < eps);
}

//// 6
TEST_CASE("Norm vector") {
    Vector A(4, 10);
    REQUIRE(abs(A.norm() - 20) < eps);
}

//// 7
TEST_CASE("Normalized vector") {
    Vector A(4, 10);
    Vector B(1);
    B = A.normalized();
    REQUIRE(B.size() == 4);
    for (int i = 0; i < B.size(); ++i) {
        REQUIRE(B[i] - 0.5 < eps);
    }
    REQUIRE(A.size() == 4);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 10) < eps);
    }
}

//// 8
TEST_CASE("Normalize vector") {
    Vector A(4, 10);
    A.normalize();
    REQUIRE(A.size() == 4);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 0.5) < eps);
    }
}

//// 9
TEST_CASE("Operator + vector") {
    Vector A(4, 10);
    Vector B(4, 5);
    Vector C(1);
    C = A + B;
    REQUIRE(C.size() == 4);
    for (int i = 0; i < C.size(); ++i) {
        REQUIRE(abs(C[i] - 15) < eps);
    }
    //////////////
    REQUIRE(A.size() == 4);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(A[i] - 10 < eps);
    }
    A += B;
    REQUIRE(A.size() == 4);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 15) < eps);
    }
}

//// 10
TEST_CASE("Operator - vector") {
    Vector A(4, 10);
    Vector B(4, 7);
    Vector C(1);
    C = A - B;
    REQUIRE(C.size() == 4);
    for (int i = 0; i < C.size(); ++i) {
        REQUIRE(abs(C[i] - 3) < eps);
    }
    //////////////
    REQUIRE(A.size() == 4);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 10) < eps);
    }
    A -= B;
    REQUIRE(A.size() == 4);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 3) < eps);
    }
}

//// 11
TEST_CASE("Operator ^ vector") {
    Vector A(4, 10);
    Vector B(4, 5);
    Vector C(1);
    C = A ^ B;
    REQUIRE(C.size() == 4);
    for (int i = 0; i < C.size(); ++i) {
        REQUIRE(abs(C[i] - 50) < eps);
    }
    //////////////
    REQUIRE(A.size() == 4);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 10) < eps);
    }
    A ^= B;
    REQUIRE(A.size() == 4);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 50) < eps);
    }
}

//// 12
TEST_CASE("Operator * vector (скалярное)") {
    Vector A(4, 10);
    Vector B(4, 5);
    REQUIRE(abs(A * B - 200) < eps);
}

//// 13
TEST_CASE("Operator * vector (vec * k)") {
    Vector A(5, 10);
    Vector B(1);
    double k = 5;
    B = A * k;
    REQUIRE(B.size() == 5);
    for (int i = 0; i < B.size(); ++i) {
        REQUIRE(abs(B[i] - 50) < eps);
    }
    ///////////////
    REQUIRE(A.size() == 5);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 10) < eps);
    }
    A *= k;
    REQUIRE(A.size() == 5);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 50) < eps);
    }
}

//// 14
TEST_CASE("Operator / vector (vec / k)") {
    Vector A(5, 10);
    Vector B(1);
    double k = 5;
    B = A / k;
    REQUIRE(B.size() == 5);
    for (int i = 0; i < B.size(); ++i) {
        REQUIRE(abs(B[i] - 2) < eps);
    }
    ///////////////
    REQUIRE(A.size() == 5);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 10) < eps);
    }
    A /= k;
    REQUIRE(A.size() == 5);
    for (int i = 0; i < A.size(); ++i) {
        REQUIRE(abs(A[i] - 2) < eps);
    }
}

//// 15
TEST_CASE("Operator == vector") {
    Vector A(5, 10);
    Vector B(5, 11);
    Vector C(4, 10);
    Vector D(5, 10);
    REQUIRE_FALSE(A == B);
    REQUIRE_FALSE(A == C);
    REQUIRE(A == D);
}

//// 16
TEST_CASE("Operator != vector") {
    Vector A(5, 10);
    Vector B(5, 11);
    Vector C(4, 10);
    Vector D(5, 10);
    REQUIRE(A != B);
    REQUIRE(A != C);
    REQUIRE_FALSE(A != D);
}

//// 17
TEST_CASE("Operator * vector (k * vec)") {
    Vector A(5, 10);
    double k = 5;
    Vector B(1);
    B = k * A;
    REQUIRE(B.size() == 5);
    for (int i = 0; i < B.size(); ++i) {
        REQUIRE(abs(B[i] - 50) < eps);
    }
}

//// 18
TEST_CASE("Operator * (vector * matrix)") {
    Matrix A(5, 3, 1);
    std::pair<size_t, size_t> p(5, 3);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            A.get(i, j) = i * p.second + j + 1;
        }
    }
    Vector B(5);
    for (int i = 0; i < B.size(); ++i) {
        B[i] = i + 1;
    }

    Vector C(1);
    C = B * A;
    REQUIRE(C.size() == 3);
    REQUIRE(abs(C[0] - 135) < eps);
    REQUIRE(abs(C[1] - 150) < eps);
    REQUIRE(abs(C[2] - 165) < eps);

    ////////////////////
    B *= A;
    REQUIRE(B.size() == 3);
    REQUIRE(abs(B[0] - 135) < eps);
    REQUIRE(abs(B[1] - 150) < eps);
    REQUIRE(abs(B[2] - 165) < eps);
}

////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MATRIX~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//// 1
TEST_CASE("Constructor matrix (square)") {
    Matrix A(5, 10.0);
    std::pair<size_t, size_t> p(5, 5);
    REQUIRE(A.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.first; ++j) {
            REQUIRE(abs(A.get(i, j) - 10) < eps);
        }
    }
}

//// 2
TEST_CASE("Eye matrix") {
    Matrix A(1);
    A = A.eye(5);
    std::pair<size_t, size_t> p(5, 5);
    REQUIRE(A.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            if (i == j) {
                REQUIRE(abs(A.get(i, j) - 1) < eps);
            }
            else {
                REQUIRE(abs(A.get(i, j)) < eps);
            }
        }
    }
}

//// 3
TEST_CASE("Constructor matrix (rows, cols)") {
    Matrix A(5, 6, 10);
    std::pair<size_t, size_t> p(5, 6);
    REQUIRE(A.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(A.get(i, j) - 10) < eps);
        }
    }
}

//// 4
TEST_CASE("Copy constructor matrix") {
    Matrix A(5, 7, 8);
    Matrix B(A);
    std::pair<size_t, size_t> p(5, 7);
    REQUIRE(B.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(B.get(i, j) - 8) < eps);
        }
    }
}

//// 5
TEST_CASE("Operator = matrix") {
    Matrix A(5, 7, 8);
    Matrix B(1);
    B = A;
    std::pair<size_t, size_t> p(5, 7);
    REQUIRE(B.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(B.get(i, j) - 8) < eps);
        }
    }
}

// 6
TEST_CASE("Reshape matrix") {
    Matrix A(5, 7, 8);
    std::pair<size_t, size_t> p(5, 7);
    REQUIRE(A.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(A.get(i, j) - 8) < eps);
        }
    }

    A.reshape(7, 5);
    std::pair<size_t, size_t> k(7, 5);
    REQUIRE(A.shape() == k);
    for (int i = 0; i < k.first; ++i) {
        for (int j = 0; j < k.second; ++j) {
            REQUIRE(abs(A.get(i, j) - 8) < eps);
        }
    }
}

//// 7
TEST_CASE("Operator [][] matrix") {
    Matrix A(5, 6, 5.0);
    double b = A.get(0, 1);
    REQUIRE(abs(b - 5) < eps);
    A.get(0, 0) = 42.0;
    b = A.get(0, 0);
    REQUIRE(abs(b - 42) < eps);
}

//// 8
TEST_CASE("Operator + matrix") {
    Matrix A(5, 7, 10.0);
    Matrix B(5, 7, 30.0);
    std::pair<size_t, size_t> p(5, 7);
    Matrix C(1);
    C = A + B;
    REQUIRE(C.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(C.get(i, j) - 40) < eps);
        }
    }
    ///////////////////
    A += B;
    REQUIRE(A.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(A.get(i, j) - 40) < eps);
        }
    }
}

//// 9
TEST_CASE("Operator - matrix") {
    Matrix A(5, 7, 30.0);
    Matrix B(5, 7, 10.0);
    std::pair<size_t, size_t> p(5, 7);
    Matrix C(1);
    C = A - B;
    REQUIRE(C.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(C.get(i, j) - 20) < eps);
        }
    }
    ///////////////////
    A -= B;
    REQUIRE(A.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(A.get(i, j) - 20) < eps);
        }
    }
}

//// 10
TEST_CASE("Operator * (matrix * matrix)") {
    Matrix A(2, 7, 5.0);
    Matrix B(7, 3, 8.0);
    std::pair<size_t, size_t> p(2, 3);
    Matrix C(1);
    C = A * B;
    REQUIRE(C.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(C.get(i, j) - 280) < eps);
        }
    }
    //////////////////
    A *= B;
    REQUIRE(A.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(A.get(i, j) - 280) < eps);
        }
    }
}

//// 11
TEST_CASE("Operator * (matrix * const)") {
    Matrix A(5, 7, 30.0);
    double k = 5;
    std::pair<size_t, size_t> p(5, 7);
    Matrix C(1);
    C = A * k;
    REQUIRE(C.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(C.get(i, j) - 150) < eps);
        }
    }
    ///////////////////
    A *= k;
    REQUIRE(A.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(A.get(i, j) - 150) < eps);
        }
    }
}

//// 12
TEST_CASE("Operator / (matrix / const") {
    Matrix A(5, 7, 30.0);
    double k = 5;
    std::pair<size_t, size_t> p(5, 7);
    Matrix C(1);
    C = A / k;
    REQUIRE(C.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(C.get(i, j) - 6) < eps);
        }
    }
    ///////////////////
    A /= k;
    REQUIRE(A.shape() == p);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            REQUIRE(abs(A.get(i, j) - 6) < eps);
        }
    }
}

//// 13
TEST_CASE("Transposed matrix, transpose matrix") {
    Matrix A(2, 3, 0);
    Matrix B(3, 2, 10);
    Matrix C(1);
    std::pair<size_t, size_t> p(2, 3);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            A.get(i, j) = i * p.second + j + 1;
        }
    }

    B.get(0, 0) = 1;
    B.get(0, 1) = 4;
    B.get(1, 0) = 2;
    B.get(1, 1) = 5;
    B.get(2, 0) = 3;
    B.get(2, 1) = 6;

    C = A.transposed();

    REQUIRE(B == C);

    A.transpose();

    REQUIRE(A == C);
}

//// 14
TEST_CASE("determinant matrix") {
    Matrix A(3, 1.0);
    A.get(0, 0) = 1;
    A.get(0, 1) = 2;
    A.get(0, 2) = 3;
    A.get(1, 0) = 4;
    A.get(1, 1) = 5;
    A.get(1, 2) = 6;
    A.get(2, 0) = 7;
    A.get(2, 1) = 8;
    A.get(2, 2) = 10;
    REQUIRE(abs(A.det() + 3) < eps);

    Matrix B(1, 5.0);
    REQUIRE(abs(B.det() - 5) < eps);
}

//// 15
TEST_CASE("Operator * (matrix * vector)") {
    Matrix A(4, 3, 1.0);
    std::pair<size_t, size_t> p(4, 3);
    for (int i = 0; i < p.first; ++i) {
        for (int j = 0; j < p.second; ++j) {
            A.get(i, j) = i * p.second + j + 1;
        }
    }

    Vector B(3);
    B[0] = 1;
    B[1] = 2;
    B[2] = 3;

    Vector C(1);
    C = A * B;

    REQUIRE(C.size() == 4);
    REQUIRE(abs(C[0] - 14) < eps);
    REQUIRE(abs(C[1] - 32) < eps);
    REQUIRE(abs(C[2] - 50) < eps);
    REQUIRE(abs(C[3] - 68) < eps);
}

//// 16
TEST_CASE("Inversion matrix") {
    Matrix A(4, 1.0);
    A.get(0, 0) = 0;
    A.get(0, 1) = 3;
    A.get(0, 2) = -1;
    A.get(0, 3) = 2;
    A.get(1, 0) = 2;
    A.get(1, 1) = 1;
    A.get(1, 2) = 0;
    A.get(1, 3) = 0;
    A.get(2, 0) = -2;
    A.get(2, 1) = -1;
    A.get(2, 2) = 0;
    A.get(2, 3) = 2;
    A.get(3, 0) = -5;
    A.get(3, 1) = 7;
    A.get(3, 2) = 1;
    A.get(3, 3) = 1;

    Matrix B(1);
    B = A.inv();

    REQUIRE(abs(B.get(0, 0) + 0.04) < eps);
    REQUIRE(abs(B.get(0, 1) - 0.46) < eps);
    REQUIRE(abs(B.get(0, 2) - 0.06) < eps);
    REQUIRE(abs(B.get(0, 3) + 0.04) < eps);
    REQUIRE(abs(B.get(1, 0) - 0.08) < eps);
    REQUIRE(abs(B.get(1, 1) - 0.08) < eps);
    REQUIRE(abs(B.get(1, 2) + 0.12) < eps);
    REQUIRE(abs(B.get(1, 3) - 0.08) < eps);
    REQUIRE(abs(B.get(2, 0) + 0.76) < eps);
    REQUIRE(abs(B.get(2, 1) - 1.24) < eps);
    REQUIRE(abs(B.get(2, 2) - 0.64) < eps);
    REQUIRE(abs(B.get(2, 3) - 0.24) < eps);
    REQUIRE(abs(B.get(3, 0)) < eps);
    REQUIRE(abs(B.get(3, 1) - 0.5) < eps);
    REQUIRE(abs(B.get(3, 2) - 0.5) < eps);
    REQUIRE(abs(B.get(3, 3)) < eps);
}

//// 17
TEST_CASE("Operator == matrix") {
    Matrix A(5, 10.0);
    Matrix B(5, 11.0);
    Matrix C(4, 10.0);
    Matrix D(5, 10.0);
    REQUIRE_FALSE(A == B);
    REQUIRE_FALSE(A == C);
    REQUIRE(A == D);
}

//// 18
TEST_CASE("Operator != matrix") {
    Matrix A(5, 10.0);
    Matrix B(5, 11.0);
    Matrix C(4, 10.0);
    Matrix D(5, 10.0);
    REQUIRE(A != B);
    REQUIRE(A != C);
    REQUIRE_FALSE(A != D);
}


