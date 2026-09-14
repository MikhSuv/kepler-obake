# AGENTS.md

## Project Overview

This project implements a C++ library of Keplerian motion functions represented as **Poisson series** (trigonometric series with polynomial coefficients). The library is intended for use in perturbation theory and celestial mechanics.

The first function to implement is `z1`. Many similar functions (`z2`, `z3`, ...) will follow. Each function will be expressed as a Poisson series in specific complex variables.

The resulting objects must support:
- Serialization to/from files (binary and/or text).
- Numerical evaluation by substituting values for variables.
- Substitution of one series into another (provided by obake).
- Truncation of the series to a given total degree in the variables (provided by obake).
- Adjustable truncation degree (stored as a parameter of the object).

## Technologies & Dependencies

- **Language**: C++20.
- **Build system**: CMake, with dependencies managed via **conda**.
- **Symbolic library**: [obake](https://github.com/bluescarni/obake) – installed via conda, linked through CMake (`find_package` or target).
- **Multiprecision**: [mp++](https://github.com/bluescarni/mppp) – rational arithmetic. Coefficients are complex numbers with rational real and imaginary parts: `std::complex<mppp::rational<1>>`. Higher precision can be obtained by choosing `mppp::rational<n>` (n > 0, fixed) or `mppp::rational<0>` (arbitrary precision).
- **Testing**: CTest (detailed test plan to be defined later).

All required dependencies are installed via a conda environment. The CMake configuration must be informed of the conda environment's prefix so that header files and libraries can be located. Typically this is done by passing the flag `-DCMAKE_PREFIX_PATH=$CONDA_PREFIX` to CMake, which automatically adds the conda include and library directories to the search paths. In some configurations it may also be necessary to add `-DCMAKE_INCLUDE_PATH=$CONDA_PREFIX/include`.

## Mathematical Background

### Keplerian variables

The basic Keplerian elements: eccentricity `e`, mean anomaly `M`, argument of pericenter `g`, longitude of ascending node `Ω`.

Define:
- `\tilde{g} = g + Ω`
- `λ = M + \tilde{g}`

Introduce complex variables and their conjugates:
- `X = sqrt(2) * sqrt(1 - sqrt(1 - e^2)) * exp(i * \tilde{g})`
- `\bar{X} = sqrt(2) * sqrt(1 - sqrt(1 - e^2)) * exp(-i * \tilde{g})`
- `Λ = exp(i * λ)`
- `\bar{Λ} = Λ^{-1} = exp(-i * λ)`

Here `i` is the imaginary unit.

In code, `Λ` is represented as a single symbolic variable (e.g. `L`); `\bar{Λ} = Λ^{-1}` is obtained as `obake::pow(L, -1)`.

### Function z1

The exact formula in Keplerian elements:
`z_1 = e*sin(M)`
In the complex variables it becomes:
`z_1 = (i/2) * [ X * Λ^{-1} - \bar{X} * Λ ] * sqrt(1 - X * \bar{X} / 4)`

To expand as a Poisson series:
1. Expand the square root as a binomial series in the variable `tmp = X * \bar{X} / 4`.
   - Use an efficient computation of binomial coefficients.
2. Multiply the expanded root by the remaining factors `(i/2) * [ X * Λ^{-1} - \bar{X} * Λ ]`.
3. Truncate the resulting series to a given degree in `X` and `\bar{X}` (the variables `Λ` and `\bar{Λ}` are treated as **parameters** that will later be substituted numerically).

*Note:* The `obake` and `mp++` libraries are not part of this repository. They are installed separately in the conda environment and are located via CMake as described above.

## Coding Conventions

- **Follow the style of the obake source code** (naming, formatting, comment style). The agent should inspect obake headers (e.g., in the conda environment) and mimic that style.
- Use C++20 features.
- All public functions and classes must be **fully documented** (Doxygen-style comments).
- Use type aliases for readability, e.g.:
  ```cpp
  using rat_t = mppp::rational<1>;   // rational numbers
  using cf_t = std::complex<rat_t>;  // coefficients
  using mono_t = obake::packed_monomial<std::int64_t>;
  using series_t = obake::p_series<mono_t, cf_t>;
  var_t = std::string; // or a custom enum for variable identification
  ```
- Ensure that all series operations (multiplication, truncation) are performed using obake's built‑in capabilities.

## Current Task: Implementation of z1
### Objective
Create a suitable symbolic object type (using obake templates) that can hold the Poisson series expansion of `z1`
with complex coefficients having rational real and imaginary parts. Show how to use this type in practice (construction, evaluation, serialization).

### Requirements for the object

- Series type: `obake::p_series<obake::packed_monomial<std::int64_t>, cf_t>` (a power series with multivariate monomial keys).
- Coefficient type: `cf_t = std::complex<mppp::rational<1>>` (or an alias) — complex numbers with rational real and imaginary parts.
- Variables: at least `X` and `\bar{X}` must be handled as symbolic variables. `Λ` is represented as a symbolic variable (e.g. `L`); `\bar{Λ}` is obtained as `Λ^{-1}` via `obake::pow(L, -1)`. `Λ` and `\bar{Λ}` are treated as parameters that will be substituted later.
- The object must support:
  - Addition, multiplication by scalar/coefficient, and multiplication of series (for the binomial expansion).
  - Substitution of a series into another (to replace `tmp = X * \bar{X} / 4 ` inside the square root expansion).
  - Truncation to a given total degree in `X` and `\bar{X}`.
  - Serialization to/from a file (binary or text).
  - Numerical evaluation given numeric values for `X`, `\bar{X}`, `Λ`, `\bar{Λ}`.

### Algorithm steps (to be implemented)
1. Construct the symbolic expression for `tmp = X * \bar{X} / 4`.
2. Expand `sqrt(1 - tmp)` as a binomial series up to a truncation degree `N` in `tmp`.
  - Use a precomputed table of binomial coefficients for efficiency.
3. Multiply the expansion by `(i/2) * [ X * Λ^{-1} - \bar{X} * Λ ]`.
4. Truncate the resulting series to degree `N` in `X` and `\bar{X}` combined (since `tmp` itself is degree 2 in these variables).
5. Store the result in the chosen obake container, together with the truncation degree and any other metadata.

### Testing plan (initial)
- Coefficient check: Verify that the first few terms of the expansion match known formulas (e.g., for small eccentricity).

- Numerical check: For a fixed eccentricity `e`, compute `z1` using the exact formula `e * sin(M)` and compare with the value obtained by evaluating the series after substituting `X`, `\bar{X}`, `Λ`, `\bar{Λ}`. Perform this for several values of `M`.

- Tests will be integrated with CTest later.
 
## Notes for the Agent
- The agent's main job is to design and implement the appropriate symbolic type using obake, and demonstrate its usage. The user will review the mathematical correctness and later integrate the object into a larger framework.

- Do not modify obake or mp++ sources; they are managed by conda and installed externally.

- If obake provides multiple candidate types (polynomial, Poisson series, etc.), the agent must compare them and justify the choice based on the requirements above.

- Keep the code modular so that the same pattern can be reused for z2, z3, etc.
