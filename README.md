# Position Based Dynamics (PBD)

A modular, high-performance C++17 physics library implementing Position Based Dynamics (PBD), Extended PBD (XPBD), and Rigid Body PBD (RPBD) algorithms for deformable bodies, rigid bodies, rods, and fluids.

## Features

- **PBD & XPBD:** Distance, isometric bending, shape matching, tetrahedral volume, and strain constraints with compliance formulations.
- **Rigid Body Dynamics (RPBD):** Ball, hinge, slider, universal, motor, and contact constraints with Delassus operator kinematics.
- **Rods:** Cosserat rods, elastic rods (discrete Darboux vectors), and stiff rods.
- **Position Based Fluids (PBF):** Incompressible fluid simulation using `CubicSpline` smoothing kernels, density constraints, and Akinci boundary handling.
- **Integrators:** Explicit Euler, symplectic Euler, and Verlet time integration schemes.
- **Dependency:** Powered by [`LinearMath`](https://github.com/Giannossk/LinearMath).

## Quick Start

### CMake Integration

```cmake
include(FetchContent)

FetchContent_Declare(
    pbd
    GIT_REPOSITORY https://github.com/Giannossk/pbd.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(pbd)

target_link_libraries(your_target PRIVATE pbd::pbd)
```

### Building & Testing

```bash
cmake -B build -DPBD_BUILD_TESTING=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

## License

Licensed under the [Apache License, Version 2.0](LICENSE).  
Copyright 2026 IOANNIS SIOKOS.
