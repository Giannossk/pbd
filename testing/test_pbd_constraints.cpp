#include "test_macros.h"
#include "pbd.h"
#include <vector>

void test_pbd_distance_constraint_equal_mass() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(2.0, 0.0, 0.0);
    Real w0 = 1.0;
    Real w1 = 1.0;
    Real restLength = 1.0;
    Real stiffness = 1.0;

    Vector3r corr0, corr1;
    bool success = pbd::solve_DistanceConstraint(p0, w0, p1, w1, restLength, stiffness, corr0, corr1);
    TEST_ASSERT(success);

    // Initial distance = 2.0, rest = 1.0, difference = +1.0
    // w0 = w1 = 1 -> delta_p0 = +0.5 along X, delta_p1 = -0.5 along X
    TEST_VEC3_NEAR(corr0, Vector3r(0.5, 0.0, 0.0), 1e-5);
    TEST_VEC3_NEAR(corr1, Vector3r(-0.5, 0.0, 0.0), 1e-5);

    Vector3r p0_new = p0 + corr0;
    Vector3r p1_new = p1 + corr1;
    Real new_length = (p0_new - p1_new).norm();
    TEST_NEAR(new_length, restLength, 1e-5);

    // Momentum conservation: m0 * corr0 + m1 * corr1 == 0
    Vector3r momentum = (1.0 / w0) * corr0 + (1.0 / w1) * corr1;
    TEST_VEC3_NEAR(momentum, Vector3r::Zero(), 1e-5);
}

void test_pbd_distance_constraint_static_particle() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(2.0, 0.0, 0.0);
    Real w0 = 0.0; // static / infinite mass
    Real w1 = 1.0; // dynamic
    Real restLength = 1.0;
    Real stiffness = 1.0;

    Vector3r corr0, corr1;
    bool success = pbd::solve_DistanceConstraint(p0, w0, p1, w1, restLength, stiffness, corr0, corr1);
    TEST_ASSERT(success);

    // Static particle should not move
    TEST_VEC3_NEAR(corr0, Vector3r::Zero(), 1e-6);
    // Dynamic particle should take all the correction (-1.0 along X)
    TEST_VEC3_NEAR(corr1, Vector3r(-1.0, 0.0, 0.0), 1e-5);

    Vector3r p1_new = p1 + corr1;
    Real new_length = (p0 - p1_new).norm();
    TEST_NEAR(new_length, restLength, 1e-5);
}

void test_pbd_distance_constraint_compression() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(0.5, 0.0, 0.0); // compressed
    Real w0 = 1.0;
    Real w1 = 1.0;
    Real restLength = 1.0;
    Real stiffness = 1.0;

    Vector3r corr0, corr1;
    bool success = pbd::solve_DistanceConstraint(p0, w0, p1, w1, restLength, stiffness, corr0, corr1);
    TEST_ASSERT(success);

    // Particles should be pushed apart
    Vector3r p0_new = p0 + corr0;
    Vector3r p1_new = p1 + corr1;
    Real new_length = (p0_new - p1_new).norm();
    TEST_NEAR(new_length, restLength, 1e-5);
}

void test_pbd_dihedral_constraint() {
    // Two triangles sharing an edge (p2, p3):
    // Triangle 0: (p0, p2, p3)
    // Triangle 1: (p1, p3, p2)
    Vector3r p2(0.0, 0.0, 0.0);
    Vector3r p3(1.0, 0.0, 0.0);
    Vector3r p0(0.5, 1.0, 0.0);
    Vector3r p1(0.5, -1.0, 0.0); // Flat in XY plane: dihedral angle is pi

    Real restAngle = 0.0;
    Real w0 = 1.0, w1 = 1.0, w2 = 1.0, w3 = 1.0;
    Real stiffness = 1.0;

    // In flat configuration, corrections should be approximately zero
    Vector3r c0, c1, c2, c3;
    bool success = pbd::solve_DihedralConstraint(p0, w0, p1, w1, p2, w2, p3, w3, restAngle, stiffness, c0, c1, c2, c3);
    TEST_ASSERT(success);
    TEST_NEAR(c0.norm(), 0.0, 1e-4);
    TEST_NEAR(c1.norm(), 0.0, 1e-4);

    // Perturb p0 out of plane
    p0.z() += static_cast<Real>(0.5);
    success = pbd::solve_DihedralConstraint(p0, w0, p1, w1, p2, w2, p3, w3, restAngle, stiffness, c0, c1, c2, c3);
    TEST_ASSERT(success);

    // Corrections should push particles to restore planarity
    TEST_ASSERT(c0.norm() > 0.0);
}

void test_pbd_volume_constraint() {
    // Tetrahedron with vertices:
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(1.0, 0.0, 0.0);
    Vector3r p2(0.0, 1.0, 0.0);
    Vector3r p3(0.0, 0.0, 1.0);

    // Volume of this right tetrahedron = 1/6
    Real restVolume = static_cast<Real>(1.0 / 6.0);
    Real w0 = 1.0, w1 = 1.0, w2 = 1.0, w3 = 1.0;
    Real stiffness = 1.0;

    // Unperturbed: volume equals rest volume, corrections ~ 0
    Vector3r c0, c1, c2, c3;
    bool success = pbd::solve_VolumeConstraint(p0, w0, p1, w1, p2, w2, p3, w3, restVolume, stiffness, c0, c1, c2, c3);
    TEST_ASSERT(success);
    TEST_NEAR(c0.norm(), 0.0, 1e-5);
    TEST_NEAR(c1.norm(), 0.0, 1e-5);
    TEST_NEAR(c2.norm(), 0.0, 1e-5);
    TEST_NEAR(c3.norm(), 0.0, 1e-5);

    // Expand tetrahedron: move p3 to (0, 0, 1.1)
    Vector3r np0 = p0, np1 = p1, np2 = p2;
    Vector3r np3(0.0, 0.0, static_cast<Real>(1.1));

    // Iterative PBD solver steps
    for (int iter = 0; iter < 15; ++iter) {
        success = pbd::solve_VolumeConstraint(np0, w0, np1, w1, np2, w2, np3, w3, restVolume, stiffness, c0, c1, c2, c3);
        TEST_ASSERT(success);
        np0 += c0;
        np1 += c1;
        np2 += c2;
        np3 += c3;
    }

    Real newVolume = (static_cast<Real>(1.0) / static_cast<Real>(6.0)) * ((np1 - np0).cross(np2 - np0)).dot(np3 - np0);
    TEST_NEAR(newVolume, restVolume, 2e-3);
}

void test_pbd_edge_point_distance_constraint() {
    Vector3r p(0.5, static_cast<Real>(0.2), 0.0);
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(1.0, 0.0, 0.0);
    Real restDist = static_cast<Real>(0.5);
    Real invMass = 1.0, invMass0 = 0.0, invMass1 = 0.0; // edge is fixed
    Real compStiffness = 1.0, stretchStiffness = 1.0;

    Vector3r corr, corr0, corr1;
    bool success = pbd::solve_EdgePointDistanceConstraint(p, invMass, p0, invMass0, p1, invMass1,
                                                          restDist, compStiffness, stretchStiffness,
                                                          corr, corr0, corr1);
    TEST_ASSERT(success);
    TEST_ASSERT(corr.norm() > 0.0);

    // Apply correction
    Vector3r p_new = p + corr;
    // Solve again: second correction should be much smaller (converging)
    Vector3r corr_next, corr0_next, corr1_next;
    pbd::solve_EdgePointDistanceConstraint(p_new, invMass, p0, invMass0, p1, invMass1,
                                          restDist, compStiffness, stretchStiffness,
                                          corr_next, corr0_next, corr1_next);
    TEST_ASSERT(corr_next.norm() < corr.norm());
}

void test_pbd_isometric_bending() {
    Vector3r p0(-1.0, 0.0, 0.0);
    Vector3r p1(1.0, 0.0, 0.0);
    Vector3r p2(0.0, 1.0, 0.0);
    Vector3r p3(0.0, -1.0, 0.0);

    Matrix4r Q;
    bool init_ok = pbd::init_IsometricBendingConstraint(p0, p1, p2, p3, Q);
    TEST_ASSERT(init_ok);

    // Perturb p0 out of the XY plane
    Vector3r p0_pert = p0 + Vector3r(0.0, 0.0, static_cast<Real>(0.4));
    Real invMass = 1.0;
    Vector3r c0, c1, c2, c3;
    bool solve_ok = pbd::solve_IsometricBendingConstraint(p0_pert, invMass, p1, invMass, p2, invMass, p3, invMass,
                                                          Q, 1.0, c0, c1, c2, c3);
    TEST_ASSERT(solve_ok);

    // Correction on p0 should push back towards z = 0
    TEST_ASSERT(c0.z() < 0.0);
}

void test_pbd_shape_matching() {
    // 8 vertices of a unit cube
    std::vector<Vector3r> x0 = {
        Vector3r(-0.5, -0.5, -0.5), Vector3r(0.5, -0.5, -0.5),
        Vector3r(0.5,  0.5, -0.5), Vector3r(-0.5, 0.5, -0.5),
        Vector3r(-0.5, -0.5,  0.5), Vector3r(0.5, -0.5,  0.5),
        Vector3r(0.5,  0.5,  0.5), Vector3r(-0.5, 0.5,  0.5)
    };
    int n = 8;
    std::vector<Real> invMasses(n, 1.0);

    Vector3r restCm;
    bool init_ok = pbd::init_ShapeMatchingConstraint(x0.data(), invMasses.data(), n, restCm);
    TEST_ASSERT(init_ok);
    TEST_VEC3_NEAR(restCm, Vector3r::Zero(), 1e-5);

    // Rotate cube by 90 degrees around Z axis and translate by (2, 3, 4)
    Matrix3r R_expected = (AngleAxisr(static_cast<Real>(3.14159265358979323846 / 2.0), Vector3r(0.0, 0.0, 1.0))).toRotationMatrix();
    Vector3r t(2.0, 3.0, 4.0);

    std::vector<Vector3r> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = R_expected * x0[i] + t;
    }

    std::vector<Vector3r> corr(n);
    Matrix3r R_extracted;
    bool solve_ok = pbd::solve_ShapeMatchingConstraint(x0.data(), x.data(), invMasses.data(), n,
                                                        restCm, 1.0, false, corr.data(), &R_extracted);
    TEST_ASSERT(solve_ok);

    // Extracted rotation should match applied rotation
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            TEST_NEAR(R_extracted(r, c), R_expected(r, c), 1e-4);
        }
    }

    // Since deformation was pure rigid motion, position corrections should be zero
    for (int i = 0; i < n; ++i) {
        TEST_NEAR(corr[i].norm(), 0.0, 1e-4);
    }

    // Now perturb one vertex: stretch x[0] outwards
    x[0] += Vector3r(0.5, 0.0, 0.0);
    solve_ok = pbd::solve_ShapeMatchingConstraint(x0.data(), x.data(), invMasses.data(), n,
                                                  restCm, 1.0, false, corr.data(), &R_extracted);
    TEST_ASSERT(solve_ok);
    // Correction on vertex 0 should pull it back (-X direction)
    TEST_ASSERT(corr[0].x() < 0.0);
}

void test_pbd_strain_triangle() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(1.0, 0.0, 0.0);
    Vector3r p2(0.0, 1.0, 0.0);

    Matrix2r invRestMat;
    bool init_ok = pbd::init_StrainTriangleConstraint(p0, p1, p2, invRestMat);
    TEST_ASSERT(init_ok);

    // Stretched triangle
    Vector3r p1_pert(static_cast<Real>(1.5), 0.0, 0.0);
    Vector3r c0, c1, c2;
    Real invMass = 1.0;
    bool solve_ok = pbd::solve_StrainTriangleConstraint(p0, invMass, p1_pert, invMass, p2, invMass,
                                                         invRestMat, 1.0, 1.0, 1.0, false, false,
                                                         c0, c1, c2);
    TEST_ASSERT(solve_ok);
    // Stretch along X should cause correction shortening the edge
    TEST_ASSERT(c1.x() < 0.0);
    TEST_ASSERT(c0.x() > 0.0);
}

void test_pbd_green_strain_and_piola_stress() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(1.0, 0.0, 0.0);
    Vector3r p2(0.0, 1.0, 0.0);
    Vector3r p3(0.0, 0.0, 1.0);

    Matrix3r invRestMat;
    Real restVolume = 0.0;
    bool init_ok = pbd::init_FEMTetraConstraint(p0, p1, p2, p3, restVolume, invRestMat);
    TEST_ASSERT(init_ok);
    TEST_NEAR(restVolume, 1.0 / 6.0, 1e-5);

    // Undeformed configuration: Green strain and Piola stress should be zero, energy = 0
    Matrix3r epsilon, sigma;
    Real energy = 0.0;
    Real mu = 1000.0, lambda = 500.0;
    pbd::computeGreenStrainAndPiolaStress(p0, p1, p2, p3, invRestMat, restVolume, mu, lambda, epsilon, sigma, energy);
    TEST_NEAR(energy, 0.0, 1e-5);

    // Deformed: expand along X
    Vector3r p1_pert(static_cast<Real>(1.2), 0.0, 0.0);
    pbd::computeGreenStrainAndPiolaStress(p0, p1_pert, p2, p3, invRestMat, restVolume, mu, lambda, epsilon, sigma, energy);
    TEST_ASSERT(energy > 0.0);
    TEST_ASSERT(sigma(0, 0) > 0.0);

    // Compute gradients
    Vector3r gradC[4];
    pbd::computeGradCGreen(restVolume, invRestMat, sigma, gradC);
    // Force sum of gradients is zero (internal forces in equilibrium)
    Vector3r sumGrad = gradC[0] + gradC[1] + gradC[2] + gradC[3];
    TEST_VEC3_NEAR(sumGrad, Vector3r::Zero(), 1e-4);
}

void test_pbd_fem_tetra() {
    Vector3r p0(0.0, 0.0, 0.0);
    Vector3r p1(1.0, 0.0, 0.0);
    Vector3r p2(0.0, 1.0, 0.0);
    Vector3r p3(0.0, 0.0, 1.0);

    Matrix3r invRestMat;
    Real restVolume = 0.0;
    pbd::init_FEMTetraConstraint(p0, p1, p2, p3, restVolume, invRestMat);

    // Perturb p1 along X
    Vector3r p1_pert(static_cast<Real>(1.3), 0.0, 0.0);
    Vector3r c0, c1, c2, c3;
    Real invMass = 1.0;
    Real youngs = 1000.0, poisson = static_cast<Real>(0.3);

    bool solve_ok = pbd::solve_FEMTetraConstraint(p0, invMass, p1_pert, invMass, p2, invMass, p3, invMass,
                                                  restVolume, invRestMat, youngs, poisson, false,
                                                  c0, c1, c2, c3);
    TEST_ASSERT(solve_ok);
    // Correction on p1 should pull back along -X
    TEST_ASSERT(c1.x() < 0.0);
}

void run_all_pbd_tests() {
    std::cout << "\n=== Running Standard PBD Constraints Tests ===" << std::endl;
    RUN_TEST(test_pbd_distance_constraint_equal_mass);
    RUN_TEST(test_pbd_distance_constraint_static_particle);
    RUN_TEST(test_pbd_distance_constraint_compression);
    RUN_TEST(test_pbd_dihedral_constraint);
    RUN_TEST(test_pbd_volume_constraint);
    RUN_TEST(test_pbd_edge_point_distance_constraint);
    RUN_TEST(test_pbd_isometric_bending);
    RUN_TEST(test_pbd_shape_matching);
    RUN_TEST(test_pbd_strain_triangle);
    RUN_TEST(test_pbd_green_strain_and_piola_stress);
    RUN_TEST(test_pbd_fem_tetra);
}
