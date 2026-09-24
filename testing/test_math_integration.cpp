#include "test_macros.h"
#include "LinearMath.h"
#include "integration.h"

void test_cross_product_matrix() {
    Vector3r a(1.0, 2.0, 3.0);
    Vector3r b(4.0, 5.0, 6.0);
    
    Matrix3r a_hat;
    matrix::crossProductMatrix(a, a_hat);

    Vector3r c_mat = a_hat * b;
    Vector3r c_cross = a.cross(b);

    TEST_VEC3_NEAR(c_mat, c_cross, 1e-5);
    
    // a_hat should be skew-symmetric: a_hat^T = -a_hat
    Matrix3r skew_sum = a_hat + a_hat.transpose();
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            TEST_NEAR(skew_sum(r, c), 0.0, 1e-6);
        }
    }
}

void test_eigen_decomposition() {
    Matrix3r A;
    A(0, 0) = 2.0; A(0, 1) = 1.0; A(0, 2) = 0.0;
    A(1, 0) = 1.0; A(1, 1) = 3.0; A(1, 2) = 1.0;
    A(2, 0) = 0.0; A(2, 1) = 1.0; A(2, 2) = 2.0;

    Matrix3r V;
    Vector3r D;
    decomposition::eigenDecomposition(A, V, D);

    // Verify V is orthogonal: V * V^T = I
    Matrix3r I = V * V.transpose();
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            Real expected = static_cast<Real>((r == c) ? 1.0 : 0.0);
            TEST_NEAR(I(r, c), expected, 1e-4);
        }
    }

    // Verify A * v_i = lambda_i * v_i
    for (int i = 0; i < 3; ++i) {
        Vector3r v_i = V.col(i);
        Vector3r Av = A * v_i;
        Vector3r lambda_v = D[i] * v_i;
        TEST_VEC3_NEAR(Av, lambda_v, 1e-4);
    }
}

void test_svd_with_inversion_handling() {
    Matrix3r A;
    A(0, 0) = 1.0; A(0, 1) = 2.0; A(0, 2) = 0.0;
    A(1, 0) = 0.0; A(1, 1) = 3.0; A(1, 2) = 4.0;
    A(2, 0) = 5.0; A(2, 1) = 0.0; A(2, 2) = 6.0;

    Matrix3r U, VT;
    Vector3r sigma;
    decomposition::svdWithInversionHandling(A, sigma, U, VT);

    // Reconstruct A = U * diag(sigma) * VT
    Matrix3r sigmaMat = Matrix3r::Zero();
    sigmaMat(0, 0) = sigma[0];
    sigmaMat(1, 1) = sigma[1];
    sigmaMat(2, 2) = sigma[2];

    Matrix3r A_reconstructed = U * sigmaMat * VT;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            TEST_NEAR(A(r, c), A_reconstructed(r, c), 1e-3);
        }
    }
}

void test_polar_decomposition() {
    Quaternionr q(AngleAxisr(static_cast<Real>(0.6), Vector3r(0.577f, 0.577f, 0.577f).normalized()));
    Matrix3r R_true = q.toRotationMatrix();

    Matrix3r S_true;
    S_true(0, 0) = 2.0; S_true(0, 1) = 0.5; S_true(0, 2) = 0.1;
    S_true(1, 0) = 0.5; S_true(1, 1) = 3.0; S_true(1, 2) = 0.2;
    S_true(2, 0) = 0.1; S_true(2, 1) = 0.2; S_true(2, 2) = 1.5;

    Matrix3r A = R_true * S_true;

    Matrix3r R, U, D;
    decomposition::polarDecomposition(A, R, U, D);

    // Verify R is orthogonal: R * R^T = I
    Matrix3r I = R * R.transpose();
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            Real expected = static_cast<Real>((r == c) ? 1.0 : 0.0);
            TEST_NEAR(I(r, c), expected, 1e-3);
        }
    }
}

void test_polar_decomposition_stable() {
    Matrix3r A;
    A(0, 0) = static_cast<Real>(1e-5); A(0, 1) = 0.0; A(0, 2) = 0.0;
    A(1, 0) = 0.0; A(1, 1) = 1.0; A(1, 2) = 0.0;
    A(2, 0) = 0.0; A(2, 1) = 0.0; A(2, 2) = 2.0;

    Matrix3r R;
    decomposition::polarDecompositionStable(A, static_cast<Real>(1e-6), R);

    // Verify R is orthogonal: R * R^T = I
    Matrix3r I = R * R.transpose();
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            Real expected = static_cast<Real>((r == c) ? 1.0 : 0.0);
            TEST_NEAR(I(r, c), expected, 1e-3);
        }
    }
}

void test_extract_rotation() {
    Quaternionr q_expected(AngleAxisr(static_cast<Real>(0.4), Vector3r(0.0, 1.0, 0.0)));
    Matrix3r R = q_expected.toRotationMatrix();

    // Dilate slightly to simulate deformation
    Matrix3r A = R * static_cast<Real>(1.05);

    Quaternionr q_extracted = Quaternionr();
    decomposition::extractRotation(A, q_extracted, 20);

    TEST_QUAT_NEAR(q_extracted, q_expected, 1e-3);
}

void test_cot_theta_and_norms() {
    Vector3r v(1.0, 0.0, 0.0);
    Vector3r w(0.0, 1.0, 0.0);

    // v and w are orthogonal, angle = 90 deg -> cot(90) = 0
    Real cot = matrix::cotTheta(v, w);
    TEST_NEAR(cot, 0.0, 1e-5);

    Matrix3r M;
    M(0, 0) = 1.0; M(0, 1) = -2.0; M(0, 2) = 3.0; // row sum = 6
    M(1, 0) = 0.0; M(1, 1) =  1.0; M(1, 2) = 1.0; // row sum = 2
    M(2, 0) = 2.0; M(2, 1) = -1.0; M(2, 2) = 1.0; // row sum = 4
    // col sums: col 0 = 3, col 1 = 4, col 2 = 5 -> oneNorm = 5
    // row sums: row 0 = 6, row 1 = 2, row 2 = 4 -> infNorm = 6

    TEST_NEAR(matrix::oneNorm(M), 5.0, 1e-5);
    TEST_NEAR(matrix::infNorm(M), 6.0, 1e-5);
}

void test_symplectic_euler_linear() {
    Real dt = static_cast<Real>(0.01);
    Real mass = static_cast<Real>(1.0);
    Vector3r pos(0.0, 1.0, 0.0);
    Vector3r vel(1.0, 0.0, 0.0);
    Vector3r accel(0.0, static_cast<Real>(-9.81), 0.0);

    integration::semiImplicitEuler(dt, mass, pos, vel, accel);
    TEST_NEAR(vel.x(), 1.0, 1e-5);
    TEST_NEAR(vel.y(), -0.0981, 1e-4);
    TEST_NEAR(pos.x(), 0.01, 1e-5);
    TEST_NEAR(pos.y(), 1.0 - 0.000981, 1e-4);

    // Test velocityUpdateFirstOrder
    Vector3r newPos = pos + Vector3r(static_cast<Real>(0.1), 0.0, 0.0);
    Vector3r updatedVel;
    integration::velocityUpdateFirstOrder(dt, mass, newPos, pos, updatedVel);
    TEST_NEAR(updatedVel.x(), 10.0, 1e-4);
}

void test_symplectic_euler_angular() {
    Real dt = static_cast<Real>(0.01);
    Real mass = static_cast<Real>(1.0);
    Matrix3r inertia = Matrix3r::Identity();
    Matrix3r invInertia = Matrix3r::Identity();
    Quaternionr rotation = Quaternionr();
    Vector3r angularVelocity(0.0, 0.0, 1.0);
    Vector3r torque = Vector3r::Zero();

    integration::semiImplicitEulerRotation(dt, mass, inertia, invInertia, rotation, angularVelocity, torque);

    // Rotation after dt around Z: angle should be approx dt * 1.0 = 0.01 rad
    TEST_NEAR(rotation.z(), static_cast<Real>(std::sin(0.005)), 1e-3);
    TEST_NEAR(rotation.w(), static_cast<Real>(std::cos(0.005)), 1e-3);

    // Test angularVelocityUpdateFirstOrder
    Quaternionr q_next = rotation;
    Vector3r computedOmega;
    Quaternionr q_identity = Quaternionr();
    integration::angularVelocityUpdateFirstOrder(dt, mass, q_next, q_identity, computedOmega);
    TEST_NEAR(computedOmega.z(), 1.0, 1e-2);
}

void test_matrix_kinematics_methods() {
    Vector3r connector(1.0, 0.5, -0.2);
    Vector3r x(0.0, 0.0, 0.0);
    Real invMass = 1.5;
    Matrix3r invInertia = Matrix3r::Identity();
    Matrix3r K;
    matrix::computeMatrixK(connector, invMass, x, invInertia, K);

    // K must be symmetric
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            TEST_NEAR(K(r, c), K(c, r), 1e-6);
        }
    }

    Quaternionr q(0.7071068f, 0.0f, 0.7071068f, 0.0f);
    Quaternionr p(0.5f, 0.5f, 0.5f, 0.5f);
    Matrix4r Q = matrix::computeMatrixQ(q);
    Matrix<Real, 4, 1> p_vec(p.w(), p.x(), p.y(), p.z());
    Matrix<Real, 4, 1> qp_vec = Q * p_vec;
    Quaternionr prod = q * p;
    TEST_NEAR(qp_vec[0], prod.w(), 1e-5);
    TEST_NEAR(qp_vec[1], prod.x(), 1e-5);
    TEST_NEAR(qp_vec[2], prod.y(), 1e-5);
    TEST_NEAR(qp_vec[3], prod.z(), 1e-5);

    Matrix4r QHat = matrix::computeMatrixQHat(q);
    Matrix<Real, 4, 1> pq_vec = QHat * p_vec;
    Quaternionr prod2 = p * q;
    TEST_NEAR(pq_vec[0], prod2.w(), 1e-5);
    TEST_NEAR(pq_vec[1], prod2.x(), 1e-5);
    TEST_NEAR(pq_vec[2], prod2.y(), 1e-5);
    TEST_NEAR(pq_vec[3], prod2.z(), 1e-5);

    Matrix<Real, 4, 3, DontAlign> G = matrix::computeMatrixG(q);
    Vector3r omega(0.0, 2.0, 0.0);
    Matrix<Real, 4, 1> q_dot = G * omega;
    Quaternionr omega_quat(0.0f, 0.0f, 2.0f, 0.0f);
    Quaternionr q_times_omega = q * omega_quat;
    TEST_NEAR(q_dot[0], static_cast<Real>(0.5) * q_times_omega.w(), 1e-5);
    TEST_NEAR(q_dot[1], static_cast<Real>(0.5) * q_times_omega.x(), 1e-5);
    TEST_NEAR(q_dot[2], static_cast<Real>(0.5) * q_times_omega.y(), 1e-5);
    TEST_NEAR(q_dot[3], static_cast<Real>(0.5) * q_times_omega.z(), 1e-5);
}

void run_all_math_integration_tests() {
    std::cout << "\n=== Running Math & Integration Tests ===" << std::endl;
    RUN_TEST(test_cross_product_matrix);
    RUN_TEST(test_eigen_decomposition);
    RUN_TEST(test_svd_with_inversion_handling);
    RUN_TEST(test_polar_decomposition);
    RUN_TEST(test_polar_decomposition_stable);
    RUN_TEST(test_extract_rotation);
    RUN_TEST(test_cot_theta_and_norms);
    RUN_TEST(test_matrix_kinematics_methods);
    RUN_TEST(test_symplectic_euler_linear);
    RUN_TEST(test_symplectic_euler_angular);
}
