#include "test_macros.h"
#include "xpbd.h"
#include "pbd.h"

void test_xpbd_distance_constraint() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(2.0, 0.0, 0.0);
    Real w0 = 1.0;
    Real w1 = 1.0;
    Real restLength = 1.0;
    Real stiffness = static_cast<Real>(1e9); // rigid constraint
    Real dt = static_cast<Real>(0.01);
    Real lambda = 0.0;

    Vector3r corr0, corr1;
    // Iterate XPBD solver steps to observe convergence
    for (int iter = 0; iter < 10; ++iter) {
        bool ok = xpbd::solve_DistanceConstraint(p0, w0, p1, w1, restLength, stiffness, dt, lambda, corr0, corr1);
        TEST_ASSERT(ok);
        p0 += corr0;
        p1 += corr1;
    }

    Real final_dist = (p0 - p1).norm();
    TEST_NEAR(final_dist, restLength, 1e-4);
    TEST_ASSERT(std::abs(lambda) > 0.0); // non-zero constraint impulse
}

void test_xpbd_distance_constraint_compliance() {
    // Test that higher compliance (lower stiffness) yields smaller position correction in a single step
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(2.0, 0.0, 0.0);
    Real restLength = 1.0;
    Real dt = static_cast<Real>(0.01);

    // Stiff
    Real lambda_stiff = 0.0;
    Vector3r c0_stiff, c1_stiff;
    xpbd::solve_DistanceConstraint(p0, 1.0, p1, 1.0, restLength, static_cast<Real>(1e6), dt, lambda_stiff, c0_stiff, c1_stiff);

    // Soft
    Real lambda_soft = 0.0;
    Vector3r c0_soft, c1_soft;
    xpbd::solve_DistanceConstraint(p0, 1.0, p1, 1.0, restLength, static_cast<Real>(10.0), dt, lambda_soft, c0_soft, c1_soft);

    // Correction of stiff spring should be greater than soft spring
    TEST_ASSERT(c0_stiff.norm() > c0_soft.norm());
}

void test_xpbd_volume_constraint() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(1.0, 0.0, 0.0);
    Vector3r p2(0.0, 1.0, 0.0);
    Vector3r p3(0.0, 0.0, static_cast<Real>(1.1));

    Real restVolume = static_cast<Real>(1.0 / 6.0);
    Real stiffness = static_cast<Real>(1e9);
    Real dt = static_cast<Real>(0.01);
    Real lambda = 0.0;

    Vector3r c0, c1, c2, c3;
    for (int iter = 0; iter < 15; ++iter) {
        bool ok = xpbd::solve_VolumeConstraint(p0, 1.0, p1, 1.0, p2, 1.0, p3, 1.0,
                                               restVolume, stiffness, dt, lambda,
                                               c0, c1, c2, c3);
        TEST_ASSERT(ok);
        p0 += c0;
        p1 += c1;
        p2 += c2;
        p3 += c3;
    }

    Real finalVolume = (static_cast<Real>(1.0) / static_cast<Real>(6.0)) * ((p1 - p0).cross(p2 - p0)).dot(p3 - p0);
    TEST_NEAR(finalVolume, restVolume, 2e-3);
}

void test_xpbd_isometric_bending() {
    Vector3r p0(-1.0, 0.0, 0.0);
    Vector3r p1(1.0, 0.0, 0.0);
    Vector3r p2(0.0, 1.0, 0.0);
    Vector3r p3(0.0, -1.0, 0.0);

    Matrix4r Q;
    bool init_ok = pbd::init_IsometricBendingConstraint(p0, p1, p2, p3, Q);
    TEST_ASSERT(init_ok);

    Vector3r p0_pert = p0 + Vector3r(0.0, 0.0, static_cast<Real>(0.3));
    Real lambda = 0.0;
    Vector3r c0, c1, c2, c3;
    bool solve_ok = xpbd::solve_IsometricBendingConstraint(p0_pert, 1.0, p1, 1.0, p2, 1.0, p3, 1.0,
                                                           Q, 1000.0, static_cast<Real>(0.01), lambda,
                                                           c0, c1, c2, c3);
    TEST_ASSERT(solve_ok);
    // Correction on p0 should push back towards flat (z < 0)
    TEST_ASSERT(c0.z() < 0.0);
}

void test_xpbd_fem_tetra() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(1.0, 0.0, 0.0);
    Vector3r p2(0.0, 1.0, 0.0);
    Vector3r p3(0.0, 0.0, 1.0);

    Matrix3r invRestMat;
    Real restVolume = 0.0;
    pbd::init_FEMTetraConstraint(p0, p1, p2, p3, restVolume, invRestMat);

    Vector3r p1_pert(static_cast<Real>(1.2), 0.0, 0.0);
    Real lambda = 0.0;
    Vector3r c0, c1, c2, c3;
    bool solve_ok = xpbd::solve_FEMTetraConstraint(p0, 1.0, p1_pert, 1.0, p2, 1.0, p3, 1.0,
                                                   restVolume, invRestMat, 10000.0, static_cast<Real>(0.3), false,
                                                   static_cast<Real>(0.01), lambda, c0, c1, c2, c3);
    TEST_ASSERT(solve_ok);
    TEST_ASSERT(c1.x() < 0.0);
}

void run_all_xpbd_tests() {
    std::cout << "\n=== Running XPBD Constraints Tests ===" << std::endl;
    RUN_TEST(test_xpbd_distance_constraint);
    RUN_TEST(test_xpbd_distance_constraint_compliance);
    RUN_TEST(test_xpbd_volume_constraint);
    RUN_TEST(test_xpbd_isometric_bending);
    RUN_TEST(test_xpbd_fem_tetra);
}
