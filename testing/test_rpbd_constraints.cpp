#include "test_macros.h"
#include "rpbd.h"

void test_rpbd_ball_joint() {
    Vector3r x0(0.0, 0.0, 0.0);
    Quaternionr q0 = Quaternionr();
    Real invMass0 = 1.0;
    Matrix3r inertia0 = Matrix3r::Identity();

    Vector3r x1(2.0, 0.0, 0.0);
    Quaternionr q1 = Quaternionr();
    Real invMass1 = 1.0;
    Matrix3r inertia1 = Matrix3r::Identity();

    Vector3r jointPos(1.0, 0.0, 0.0);

    Matrix<Real, 3, 4, DontAlign> jointInfo;
    bool init_ok = rpbd::init_BallJoint(x0, q0, x1, q1, jointPos, jointInfo);
    TEST_ASSERT(init_ok);

    // Perturb body 1 position: separate joint
    x1 += Vector3r(static_cast<Real>(0.3), static_cast<Real>(0.2), static_cast<Real>(-0.1));

    Vector3r corr_x0, corr_x1;
    Quaternionr corr_q0, corr_q1;

    for (int iter = 0; iter < 5; ++iter) {
        rpbd::update_BallJoint(x0, q0, x1, q1, jointInfo);
        bool solve_ok = rpbd::solve_BallJoint(invMass0, x0, inertia0, q0,
                                              invMass1, x1, inertia1, q1,
                                              jointInfo, corr_x0, corr_q0, corr_x1, corr_q1);
        TEST_ASSERT(solve_ok);

        // Additive position and quaternion updates in PBD
        x0 += corr_x0;
        x1 += corr_x1;
        q0.coeffs() += corr_q0.coeffs();
        q0.normalize();
        q1.coeffs() += corr_q1.coeffs();
        q1.normalize();
    }

    rpbd::update_BallJoint(x0, q0, x1, q1, jointInfo);
    Vector3r c0_world = jointInfo.col(2);
    Vector3r c1_world = jointInfo.col(3);

    // Connectors should coincide after solve
    TEST_VEC3_NEAR(c0_world, c1_world, 1e-4);
}

void test_rpbd_ball_on_line_joint() {
    Vector3r x0(0.0, 0.0, 0.0);
    Quaternionr q0 = Quaternionr();
    Vector3r x1(2.0, 0.0, 0.0);
    Quaternionr q1 = Quaternionr();
    Vector3r pos(1.0, 0.0, 0.0);
    Vector3r dir(1.0, 0.0, 0.0); // Line along X axis

    Matrix<Real, 3, 10, DontAlign> jointInfo;
    bool init_ok = rpbd::init_BallOnLineJoint(x0, q0, x1, q1, pos, dir, jointInfo);
    TEST_ASSERT(init_ok);

    // Displace body 1 in Y (perpendicular to line) and along X (along line)
    x1 += Vector3r(static_cast<Real>(0.5), static_cast<Real>(0.3), 0.0);

    Vector3r corr_x0, corr_x1;
    Quaternionr corr_q0, corr_q1;
    Matrix3r I = Matrix3r::Identity();

    for (int iter = 0; iter < 5; ++iter) {
        rpbd::update_BallOnLineJoint(x0, q0, x1, q1, jointInfo);
        bool solve_ok = rpbd::solve_BallOnLineJoint(1.0, x0, I, q0, 1.0, x1, I, q1,
                                                    jointInfo, corr_x0, corr_q0, corr_x1, corr_q1);
        TEST_ASSERT(solve_ok);

        x0 += corr_x0;
        x1 += corr_x1;
        q0.coeffs() += corr_q0.coeffs();
        q0.normalize();
        q1.coeffs() += corr_q1.coeffs();
        q1.normalize();
    }

    rpbd::update_BallOnLineJoint(x0, q0, x1, q1, jointInfo);
    Vector3r c0 = jointInfo.col(5);
    Vector3r c1 = jointInfo.col(6);
    Vector3r diff = c1 - c0;

    // Perpendicular component (Y and Z) should be resolved
    TEST_NEAR(diff.y(), 0.0, 1e-4);
    TEST_NEAR(diff.z(), 0.0, 1e-4);
}

void test_rpbd_hinge_joint() {
    Vector3r x0(0.0, 0.0, 0.0);
    Quaternionr q0 = Quaternionr();
    Vector3r x1(2.0, 0.0, 0.0);
    Quaternionr q1 = Quaternionr();
    Vector3r pos(1.0, 0.0, 0.0);
    Vector3r axis(0.0, 0.0, 1.0); // Hinge around Z axis

    Matrix<Real, 4, 7, DontAlign> jointInfo;
    bool init_ok = rpbd::init_HingeJoint(x0, q0, x1, q1, pos, axis, jointInfo);
    TEST_ASSERT(init_ok);

    // Displace body 1 and rotate it off-axis (around X axis)
    x1 += Vector3r(static_cast<Real>(0.1), static_cast<Real>(0.2), 0.0);
    q1 = Quaternionr(AngleAxisr(static_cast<Real>(0.1), Vector3r(1.0, 0.0, 0.0)));

    Vector3r corr_x0, corr_x1;
    Quaternionr corr_q0, corr_q1;
    Matrix3r I = Matrix3r::Identity();

    for (int iter = 0; iter < 5; ++iter) {
        rpbd::update_HingeJoint(x0, q0, x1, q1, jointInfo);
        bool solve_ok = rpbd::solve_HingeJoint(1.0, x0, I, q0, 1.0, x1, I, q1,
                                               jointInfo, corr_x0, corr_q0, corr_x1, corr_q1);
        TEST_ASSERT(solve_ok);

        x0 += corr_x0;
        x1 += corr_x1;
        q0.coeffs() += corr_q0.coeffs();
        q0.normalize();
        q1.coeffs() += corr_q1.coeffs();
        q1.normalize();
    }

    rpbd::update_HingeJoint(x0, q0, x1, q1, jointInfo);
    Vector3r c0(jointInfo(0, 4), jointInfo(1, 4), jointInfo(2, 4));
    Vector3r c1(jointInfo(0, 5), jointInfo(1, 5), jointInfo(2, 5));
    TEST_VEC3_NEAR(c0, c1, 1e-4);
}

void test_rpbd_slider_joint() {
    Vector3r x0(0.0, 0.0, 0.0);
    Quaternionr q0 = Quaternionr();
    Vector3r x1(2.0, 0.0, 0.0);
    Quaternionr q1 = Quaternionr();
    Vector3r axis(1.0, 0.0, 0.0); // Slider along X axis

    Matrix<Real, 4, 6, DontAlign> jointInfo;
    bool init_ok = rpbd::init_SliderJoint(x0, q0, x1, q1, axis, jointInfo);
    TEST_ASSERT(init_ok);

    // 1. Test translation: displace body 1 along slide axis (0.5) and perpendicular (0.2 in Y)
    x1 += Vector3r(static_cast<Real>(0.5), static_cast<Real>(0.2), 0.0);

    Vector3r corr_x0, corr_x1;
    Quaternionr corr_q0, corr_q1;
    Matrix3r I = Matrix3r::Identity();

    rpbd::update_SliderJoint(x0, q0, x1, q1, jointInfo);
    bool solve_ok = rpbd::solve_SliderJoint(1.0, x0, I, q0, 1.0, x1, I, q1,
                                            jointInfo, corr_x0, corr_q0, corr_x1, corr_q1);
    TEST_ASSERT(solve_ok);

    x0 += corr_x0;
    x1 += corr_x1;

    // Perpendicular distance between centers should be zero
    Vector3r diff = x1 - x0;
    TEST_NEAR(diff.y(), 0.0, 1e-4);
    TEST_NEAR(diff.z(), 0.0, 1e-4);
    // Sliding along X axis should be preserved (initial 2.0 + 0.5 = 2.5)
    TEST_NEAR(diff.x(), static_cast<Real>(2.5), 1e-4);

    // 2. Test rotational locking: perturb orientation
    q1 = Quaternionr(AngleAxisr(static_cast<Real>(0.05), Vector3r(0.0, 0.0, 1.0)));
    for (int iter = 0; iter < 10; ++iter) {
        rpbd::update_SliderJoint(x0, q0, x1, q1, jointInfo);
        rpbd::solve_SliderJoint(1.0, x0, I, q0, 1.0, x1, I, q1,
                                jointInfo, corr_x0, corr_q0, corr_x1, corr_q1);
        q0.coeffs() += corr_q0.coeffs();
        q0.normalize();
        q1.coeffs() += corr_q1.coeffs();
        q1.normalize();
    }
    Quaternionr q_rel = q0.inverse() * q1;
    TEST_QUAT_NEAR(q_rel, Quaternionr(), 1e-3);
}

void run_all_rpbd_tests() {
    std::cout << "\n=== Running RPBD Joint Constraints Tests ===" << std::endl;
    RUN_TEST(test_rpbd_ball_joint);
    RUN_TEST(test_rpbd_ball_on_line_joint);
    RUN_TEST(test_rpbd_hinge_joint);
    RUN_TEST(test_rpbd_slider_joint);
}
