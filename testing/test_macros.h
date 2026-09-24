#pragma once

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>

#include "LinearMath.h"

inline int g_tests_passed = 0;
inline int g_tests_failed = 0;

#define TEST_ASSERT(cond) \
    do { \
        if (!(cond)) { \
            std::cerr << "\n    FAILED: " #cond " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            g_tests_failed++; \
            return; \
        } \
    } while (0)

#define TEST_NEAR(a, b, eps) \
    do { \
        LinearMath::Real diff_val = std::abs(static_cast<LinearMath::Real>(a) - static_cast<LinearMath::Real>(b)); \
        if (diff_val > static_cast<LinearMath::Real>(eps)) { \
            std::cerr << "\n    FAILED: " #a " (" << (a) << ") vs " #b " (" << (b) \
                      << ") diff=" << diff_val << " > " << (eps) \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            g_tests_failed++; \
            return; \
        } \
    } while (0)

#define TEST_VEC3_NEAR(v1, v2, eps) \
    do { \
        LinearMath::Real diff_vec = ((v1) - (v2)).norm(); \
        if (diff_vec > static_cast<LinearMath::Real>(eps)) { \
            std::cerr << "\n    FAILED: Vector3r near check: (" \
                      << (v1).x() << ", " << (v1).y() << ", " << (v1).z() << ") vs (" \
                      << (v2).x() << ", " << (v2).y() << ", " << (v2).z() << ") diff=" \
                      << diff_vec << " > " << (eps) \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            g_tests_failed++; \
            return; \
        } \
    } while (0)

#define TEST_QUAT_NEAR(q1, q2, eps) \
    do { \
        LinearMath::Real d1 = ((q1).coeffs() - (q2).coeffs()).norm(); \
        LinearMath::Real d2 = ((q1).coeffs() + (q2).coeffs()).norm(); \
        LinearMath::Real q_diff = std::min(d1, d2); \
        if (q_diff > static_cast<LinearMath::Real>(eps)) { \
            std::cerr << "\n    FAILED: Quaternion near check: diff=" \
                      << q_diff << " > " << (eps) \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            g_tests_failed++; \
            return; \
        } \
    } while (0)

#define RUN_TEST(test_func) \
    do { \
        int before = g_tests_failed; \
        std::cout << "  [TEST] " #test_func "... " << std::flush; \
        test_func(); \
        if (g_tests_failed == before) { \
            std::cout << "PASSED" << std::endl; \
            g_tests_passed++; \
        } \
    } while (0)
