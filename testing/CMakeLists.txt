add_executable(test_pbd
    main.cpp
    test_math_integration.cpp
    test_pbd_constraints.cpp
    test_xpbd_constraints.cpp
    test_rpbd_constraints.cpp
    test_rods_constraints.cpp
    test_fluids_constraints.cpp
)

target_link_libraries(test_pbd PRIVATE pbd LinearMath)
target_include_directories(test_pbd PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
)

if(MSVC)
    target_compile_options(test_pbd PRIVATE /EHsc /W4)
else()
    target_compile_options(test_pbd PRIVATE -Wall -Wextra)
endif()

add_test(NAME Math_Integration COMMAND test_pbd --math)
add_test(NAME PBD_Standard_Constraints COMMAND test_pbd --pbd)
add_test(NAME XPBD_Constraints COMMAND test_pbd --xpbd)
add_test(NAME RPBD_Joint_Constraints COMMAND test_pbd --rpbd)
add_test(NAME Rod_Constraints COMMAND test_pbd --rods)
add_test(NAME Fluid_Constraints COMMAND test_pbd --fluids)
add_test(NAME All_Constraints_Test COMMAND test_pbd)
