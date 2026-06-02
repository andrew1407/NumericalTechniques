# NumericalTechniques

A collection of standalone C++ programs implementing classic numerical methods, organized as a series of self-contained lab exercises. Each program embeds its own input data and prints the iterations and results to the console.

## Methods Implemented

| File | Method |
| --- | --- |
| `sripts/2-lab.cpp` | Seidel iteration method for systems of linear equations |
| `sripts/3-lab.cpp` | Chord (secant) method for root finding |
| `sripts/4.1-lab.cpp` | Simplified Newton method for nonlinear systems |
| `sripts/4.2-lab.cpp` | Steepest descent method for nonlinear systems |
| `sripts/5-lab.cpp` | Table function interpolation (divided and finite differences) |
| `sripts/6-lab.cpp` | Table function approximation (least squares) |
| `sripts/7-lab.cpp` | Numerical integration (trapezoidal and Simpson's rule) |
| `sripts/8-lab.cpp` | ODE solving (fourth-order Runge-Kutta and Euler methods) |

## Tech Stack

- C++ (standard library and `<math.h>`)
- A `c++` compiler (e.g. GCC or Clang)
- Bash (for the helper run script)

## Build and Run

Each program is independent. Compile and run any one directly:

```bash
c++ sripts/2-lab.cpp -o script
./script
```

Alternatively, use the `lab` helper script in the `sripts` directory, which compiles the chosen lab, runs it, and removes the temporary binary. Pass the lab number (the part before `-lab.cpp`):

```bash
cd sripts
./lab 2     # builds and runs 2-lab.cpp
```

## Project Structure

```
.
├── .gitignore
└── sripts/
    ├── lab            # bash helper: compiles, runs, and cleans up a lab
    ├── 2-lab.cpp      # Seidel iteration method
    ├── 3-lab.cpp      # chord method
    ├── 4.1-lab.cpp    # simplified Newton method
    ├── 4.2-lab.cpp    # steepest descent method
    ├── 5-lab.cpp      # table function interpolation
    ├── 6-lab.cpp      # table function approximation
    ├── 7-lab.cpp      # numerical integration
    └── 8-lab.cpp      # Runge-Kutta and Euler methods
```
