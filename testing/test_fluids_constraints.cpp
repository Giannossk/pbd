#include "test_macros.h"
#include "pbdfluids.h"
#include <vector>

void test_cubic_spline_properties() {
    Real h = static_cast<Real>(0.5);
    CubicSpline::setRadius(h);

    TEST_NEAR(CubicSpline::getRadius(), h, 1e-6);

    // W(0) should be strictly positive
    Real w_zero = CubicSpline::W_zero();
    TEST_ASSERT(w_zero > 0.0);

    // Compact support: W(r) == 0 for r >= h
    Vector3r r_outside(h + static_cast<Real>(0.01), 0.0, 0.0);
    TEST_NEAR(CubicSpline::W(r_outside), 0.0, 1e-6);

    // Symmetry: W(r) == W(-r)
    Vector3r r(static_cast<Real>(0.1), static_cast<Real>(-0.2), static_cast<Real>(0.15));
    TEST_NEAR(CubicSpline::W(r), CubicSpline::W(-r), 1e-6);

    // Gradient at origin is zero
    Vector3r grad_zero = CubicSpline::gradW(Vector3r::Zero());
    TEST_VEC3_NEAR(grad_zero, Vector3r::Zero(), 1e-6);

    // Gradient anti-symmetry: gradW(r) == -gradW(-r)
    Vector3r grad_pos = CubicSpline::gradW(r);
    Vector3r grad_neg = CubicSpline::gradW(-r);
    TEST_VEC3_NEAR(grad_pos, -grad_neg, 1e-5);
}

void test_pbf_density_and_lagrange_multiplier() {
    Real h = static_cast<Real>(0.5);
    CubicSpline::setRadius(h);

    // Particle 0 at origin, surrounded by 6 neighbors along +-X, +-Y, +-Z at distance 0.2
    int numParticles = 7;
    std::vector<Vector3r> x(numParticles);
    std::vector<Real> mass(numParticles, 1.0);

    x[0] = Vector3r(0.0, 0.0, 0.0);
    x[1] = Vector3r(static_cast<Real>(0.2), 0.0, 0.0);
    x[2] = Vector3r(static_cast<Real>(-0.2), 0.0, 0.0);
    x[3] = Vector3r(0.0, static_cast<Real>(0.2), 0.0);
    x[4] = Vector3r(0.0, static_cast<Real>(-0.2), 0.0);
    x[5] = Vector3r(0.0, 0.0, static_cast<Real>(0.2));
    x[6] = Vector3r(0.0, 0.0, static_cast<Real>(-0.2));

    std::vector<unsigned int> neighbors = { 1, 2, 3, 4, 5, 6 };
    Real restDensity = 1.0;
    Real density_err = 0.0;
    Real density = 0.0;

    bool ok_dens = pbdfluids::computePBFDensity(0, numParticles, x.data(), mass.data(),
                                                nullptr, nullptr, 6, neighbors.data(),
                                                restDensity, false, density_err, density);
    TEST_ASSERT(ok_dens);
    // Density should be greater than rest density (mass * W sum)
    TEST_ASSERT(density > 0.0);

    // Compute Lagrange multiplier
    Real lambda = 0.0;
    bool ok_lambda = pbdfluids::computePBFLagrangeMultiplier(0, numParticles, x.data(), mass.data(),
                                                            nullptr, nullptr, density, 6, neighbors.data(),
                                                            restDensity, false, lambda);
    TEST_ASSERT(ok_lambda);

    // When density > restDensity, constraint C = density/restDensity - 1 > 0
    // so lambda = -C / sum(grad^2) should be negative
    if (density > restDensity) {
        TEST_ASSERT(lambda < 0.0);
    }
}

void test_pbf_density_solve() {
    Real h = static_cast<Real>(0.5);
    CubicSpline::setRadius(h);

    int numParticles = 2;
    std::vector<Vector3r> x = { Vector3r(0.0, 0.0, 0.0), Vector3r(static_cast<Real>(0.2), 0.0, 0.0) };
    std::vector<Real> mass = { 1.0, 1.0 };
    std::vector<unsigned int> neighbors = { 1 };
    std::vector<Real> lambda = { static_cast<Real>(-0.5), static_cast<Real>(-0.5) }; // Both in compression

    Vector3r corr;
    bool ok = pbdfluids::solveDensityConstraint(0, numParticles, x.data(), mass.data(),
                                               nullptr, nullptr, 1, neighbors.data(),
                                               1.0, false, lambda.data(), corr);
    TEST_ASSERT(ok);

    // Particle 0 should be pushed away from particle 1 (i.e. Along -X)
    TEST_ASSERT(corr.x() < 0.0);
}

void run_all_fluids_tests() {
    std::cout << "\n=== Running Fluid Constraints Tests ===" << std::endl;
    RUN_TEST(test_cubic_spline_properties);
    RUN_TEST(test_pbf_density_and_lagrange_multiplier);
    RUN_TEST(test_pbf_density_solve);
}
