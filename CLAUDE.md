# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PCbound is a set of standalone C programs for modeling pulsar polar cap boundary dynamics — simulating spark configurations, drift motion, corotation effects, and line-of-sight emission. Each `.c` file compiles into an independent executable.

## Build Commands

```sh
make        # build all programs into ~/bin/
make clean  # remove all built binaries
```

Each program can also be built individually, e.g. `make ~/bin/PCdrift`.

Individual compile commands are also embedded in a comment at the top of each `.c` file. General pattern:
```sh
gcc <source>.c -o ~/bin/<name> [-lcpgplot] [-lfftw3] -lm
```

External dependencies: `libfftw3-dev` (FFT), `giza-dev` (provides `libcpgplot`), `libm`.

## Data Pipeline

Programs are meant to be run in sequence; outputs of earlier steps feed into later ones:

```
polcap infield outopen alpha
    → produces open field boundary file

tracklos infield inopen trkfile alpha beta nlos
    → produces line-of-sight track file

PCdrift inopen trkfile outfile alpha beta P3 npulse
    → main simulation: single pulses with drifting sparks + LRFS analysis
```

## Architecture

Each `.c` file is a self-contained program with a `main()`. Shared logic lives in header files that are `#include`d directly (no separate compilation units):

- **sparkfunc.h** — 3D vector algebra (`unitvect`, `crossprod`, `dotprod`, `rotvect`) and `sparkconfig()` for placing sparks on the polar cap
- **usefunc.h** — General utilities: `urandom()` (Park-Miller RNG), `meanrms()`, `noise()` (Gaussian noise), `rec2pol()`, `coortrans()` (coordinate rotation)
- **plotfunc.h** — FFT via FFTW3 (`pulsfft()`) and PGPLOT visualization (`plotlrfs()`)
- **ellipse_fit.h** — Least-squares ellipse fitting (`FitEllipse()`, `LeastSquares()`, `LinearSystemByGaussian()`)
- **gauss_elim.h** — Standalone Gaussian elimination with partial pivoting
- **PCdrift.h** — Function declarations for PCdrift.c

## No Tests or Linting

There is no test suite and no linter/formatter configuration.
