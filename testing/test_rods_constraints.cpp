#include "test_macros.h"
#include "pbdelasticrods.h"
#include "pbdcosseratrods.h"

void test_elastic_rods_perpendicular_bisector() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(2.0, 0.0, 0.0);
    // Midpoint is at (1, 0, 0).
    // Vector p1 - p0 = (2, 0, 0).
    // Perpendicular point p2 should have x = 1.0 (so p2 - pm is orthogonal to p1 - p0).
    // Let's set p2 at (1.5, 1.0, 0.0), which violates perpendicularity.
    Vector3r p2(1.5, 1.0, 0.0);

    Real invMass0 = 1.0, invMass1 = 1.0, invMass2 = 1.0;
    Real stiffness = 1.0;
    Vector3r corr0, corr1, corr2;

    Vector3r np0 = p0, np1 = p1, np2 = p2;
    for (int iter = 0; iter < 5; ++iter) {
        bool ok = pbdelasticrods::solve_PerpendiculaBisectorConstraint(np0, invMass0, np1, invMass1, np2, invMass2,
                                                                       stiffness, corr0, corr1, corr2);
        TEST_ASSERT(ok);
        np0 += corr0;
        np1 += corr1;
        np2 += corr2;
    }

    Vector3r npm = 0.5 * (np0 + np1);
    Real dotProd = (np2 - npm).dot(np1 - np0);
    TEST_NEAR(dotProd, 0.0, 1e-4);
}

void test_elastic_rods_ghost_point_distance() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(2.0, 0.0, 0.0);
    // Midpoint = (1, 0, 0).
    // Point p2 at (1.0, 0.5, 0.0) -> current distance = 0.5.
    Vector3r p2(1.0, 0.5, 0.0);
    Real restDist = 1.0;

    Vector3r c0, c1, c2;
    bool ok = pbdelasticrods::solve_GhostPointEdgeDistanceConstraint(p0, 1.0, p1, 1.0, p2, 1.0,
                                                                     1.0, restDist, c0, c1, c2);
    TEST_ASSERT(ok);

    Vector3r np0 = p0 + c0;
    Vector3r np1 = p1 + c1;
    Vector3r np2 = p2 + c2;
    Vector3r npm = 0.5 * (np0 + np1);
    Real newDist = (np2 - npm).norm();
    TEST_NEAR(newDist, restDist, 1e-4);
}

void test_elastic_rods_material_frame() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(1.0, 0.0, 0.0);
    Vector3r p2(0.5, 1.0, 0.0); // Ghost point above midpoint

    Matrix3r frame;
    bool ok = pbdelasticrods::computeMaterialFrame(p0, p1, p2, frame);
    TEST_ASSERT(ok);

    // Frame must be orthonormal: columns have unit length and are orthogonal
    for (int i = 0; i < 3; ++i) {
        TEST_NEAR(frame.col(i).norm(), 1.0, 1e-5);
    }
    TEST_NEAR(frame.col(0).dot(frame.col(1)), 0.0, 1e-5);
    TEST_NEAR(frame.col(0).dot(frame.col(2)), 0.0, 1e-5);
    TEST_NEAR(frame.col(1).dot(frame.col(2)), 0.0, 1e-5);

    // Determinant should be +1 (right-handed frame)
    TEST_NEAR(frame.determinant(), 1.0, 1e-5);
}

void test_elastic_rods_darboux_vector() {
    Vector3r p0(0.0, 0.0, 0.0), p1(1.0, 0.0, 0.0), p2(0.5, 1.0, 0.0);
    Matrix3r frameA;
    pbdelasticrods::computeMaterialFrame(p0, p1, p2, frameA);

    // Between identical frames, Darboux vector must be zero
    Vector3r darboux;
    bool ok = pbdelasticrods::computeDarbouxVector(frameA, frameA, 1.0, darboux);
    TEST_ASSERT(ok);
    TEST_VEC3_NEAR(darboux, Vector3r::Zero(), 1e-5);
}

void test_cosserat_rods_stretch_shear() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(1.5, 0.0, 0.0); // Stretched edge (rest = 1.0)
    Quaternionr q0 = Quaternionr(); // Tangent along (1, 0, 0) in local coords
    Real restLength = 1.0;
    Vector3r ks(100.0, 100.0, 100.0);

    Vector3r corr0, corr1;
    Quaternionr corrq0;
    bool ok = pbdcosseratrods::solve_StretchShearConstraint(p0, 1.0, p1, 1.0, q0, 1.0,
                                                           ks, restLength, corr0, corr1, corrq0);
    TEST_ASSERT(ok);

    // Particles should be pulled together along X
    TEST_ASSERT(corr0.x() > 0.0);
    TEST_ASSERT(corr1.x() < 0.0);
}

void test_cosserat_rods_bend_twist() {
    Quaternionr q0 = Quaternionr();
    // Rotate q1 slightly around Y
    Quaternionr q1(AngleAxisr(static_cast<Real>(0.1), Vector3r(0.0, 1.0, 0.0)));
    Quaternionr restDarboux = Quaternionr(); // Rest is parallel
    Vector3r ks(100.0, 100.0, 100.0);

    Quaternionr corrq0, corrq1;
    bool ok = pbdcosseratrods::solve_BendTwistConstraint(q0, 1.0, q1, 1.0,
                                                         ks, restDarboux, corrq0, corrq1);
    TEST_ASSERT(ok);
    // Corrections should rotate them towards each other
    TEST_ASSERT(corrq0.coeffs().norm() > 0.0 || corrq1.coeffs().norm() > 0.0);
}

void run_all_rods_tests() {
    std::cout << "\n=== Running Rod Constraints Tests ===" << std::endl;
    RUN_TEST(test_elastic_rods_perpendicular_bisector);
    RUN_TEST(test_elastic_rods_ghost_point_distance);
    RUN_TEST(test_elastic_rods_material_frame);
    RUN_TEST(test_elastic_rods_darboux_vector);
    RUN_TEST(test_cosserat_rods_stretch_shear);
    RUN_TEST(test_cosserat_rods_bend_twist);
}
