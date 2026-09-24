#include <iostream>
#include <string>
#include <vector>

#include "test_macros.h"

// Forward declarations
void run_all_math_integration_tests();
void run_all_pbd_tests();
void run_all_xpbd_tests();
void run_all_rpbd_tests();
void run_all_rods_tests();
void run_all_fluids_tests();

int main(int argc, char* argv[]) {
    bool run_math = false;
    bool run_pbd = false;
    bool run_xpbd = false;
    bool run_rpbd = false;
    bool run_rods = false;
    bool run_fluids = false;

    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--math") run_math = true;
            else if (arg == "--pbd") run_pbd = true;
            else if (arg == "--xpbd") run_xpbd = true;
            else if (arg == "--rpbd") run_rpbd = true;
            else if (arg == "--rods") run_rods = true;
            else if (arg == "--fluids") run_fluids = true;
        }
    } else {
        // Run all by default
        run_math = true;
        run_pbd = true;
        run_xpbd = true;
        run_rpbd = true;
        run_rods = true;
        run_fluids = true;
    }

    std::cout << "========================================\n"
              << "       PBD Constraint Math Test Suite   \n"
              << "========================================" << std::endl;

    if (run_math)   run_all_math_integration_tests();
    if (run_pbd)    run_all_pbd_tests();
    if (run_xpbd)   run_all_xpbd_tests();
    if (run_rpbd)   run_all_rpbd_tests();
    if (run_rods)   run_all_rods_tests();
    if (run_fluids) run_all_fluids_tests();

    std::cout << "\n========================================\n"
              << "              TEST SUMMARY              \n"
              << "========================================\n"
              << "  Total Passed: " << g_tests_passed << "\n"
              << "  Total Failed: " << g_tests_failed << "\n"
              << "========================================" << std::endl;

    if (g_tests_failed > 0) {
        std::cerr << "OVERALL STATUS: FAILURE" << std::endl;
        return 1;
    }

    std::cout << "OVERALL STATUS: ALL TESTS PASSED!" << std::endl;
    return 0;
}
