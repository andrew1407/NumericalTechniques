# NumericalTechniques

<p align="center">
  <a href="https://andrew1407.github.io/NumericalTechniques/" title="Open the labs in the browser">
    <img src="web/icon.svg" alt="Numerical Techniques Labs" width="140" height="140">
  </a>
</p>

[![pages deploy badge](https://github.com/andrew1407/NumericalTechniques/actions/workflows/pages.yml/badge.svg)](https://andrew1407.github.io/NumericalTechniques/)

A collection of standalone C++ programs implementing classic numerical methods, organized as a series of self-contained lab exercises. Each program embeds its own input data and prints the iterations and results to the console.

## Methods Implemented

| File | Method |
| --- | --- |
| [`sripts/2-lab.cpp`](sripts/2-lab.cpp) | Seidel iteration method for systems of linear equations |
| [`sripts/3-lab.cpp`](sripts/3-lab.cpp) | Chord (secant) method for root finding |
| [`sripts/4.1-lab.cpp`](sripts/4.1-lab.cpp) | Simplified Newton method for nonlinear systems |
| [`sripts/4.2-lab.cpp`](sripts/4.2-lab.cpp) | Steepest descent method for nonlinear systems |
| [`sripts/5-lab.cpp`](sripts/5-lab.cpp) | Table function interpolation (divided and finite differences) |
| [`sripts/6-lab.cpp`](sripts/6-lab.cpp) | Table function approximation (least squares) |
| [`sripts/7-lab.cpp`](sripts/7-lab.cpp) | Numerical integration (trapezoidal and Simpson's rule) |
| [`sripts/8-lab.cpp`](sripts/8-lab.cpp) | ODE solving (fourth-order Runge-Kutta and Euler methods) |

## Showcase

**[Open the labs in the browser](https://andrew1407.github.io/NumericalTechniques/)**. Every lab runs natively: each `sripts/*-lab.cpp` is compiled unchanged to WebAssembly with [Emscripten](https://emscripten.org), and the page runs its `main()` in a web worker, then turns what it prints into tables and interactive charts (zoom, log scale, series toggles). The raw console output is shown next to them, and every chart (PNG or CSV), table and output can be copied.

- The labs read no input, so their constants stay as compiled. Labs 2 and 5 keep their data in mutable globals (`M`, and `X`, `Y`), so the page lets you edit the system or the table and writes it into the program's memory before `main()` starts.
- Curves badged **JS** (the exact solution of the linear system, the root of the nonlinear system, a reference integral) are computed by the page as an independent check.

Run it locally (needs [emsdk](https://emscripten.org/docs/getting_started/downloads.html)):

```sh
bash web/build.sh                                  # builds web/wasm/
python -m http.server                              # then open http://localhost:8000/web/
```

The [pages workflow](.github/workflows/pages.yml) builds and deploys it on every push to `main`.

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

Alternatively, use the [`lab`](sripts/lab) helper script in the [`sripts`](sripts) directory, which compiles the chosen lab, runs it, and removes the temporary binary. Pass the lab number (the part before `-lab.cpp`):

```bash
cd sripts
./lab 2     # builds and runs 2-lab.cpp
```

## Project Structure

```
.
├── .gitignore
├── .github/workflows/pages.yml  # builds the WebAssembly modules and deploys the showcase
├── web/               # in-browser showcase: index.html, app.js, worker.js, build.sh (emscripten)
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
