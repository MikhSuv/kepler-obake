// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef KEPLER_EVALUATE_HPP
#define KEPLER_EVALUATE_HPP

#include <cstddef>
#include <vector>

#include <mp++/complex.hpp>
#include <mp++/real.hpp>

#include <obake/symbols.hpp>

#include "series.hpp"
#include "types.hpp"

namespace kepler
{

/// Numerically evaluate a Keplerian series at a given set of numeric values
/// for its symbols.
///
/// Wraps the per-term evaluation of the underlying obake series, converting
/// the complex rational coefficients to the arbitrary-precision complex type
/// \ref cnum_t.
[[nodiscard]] cnum_t evaluate(const series &s, const obake::symbol_map<cnum_t> &sm);

/// Evaluate a Keplerian series on a uniform grid of the angle \f$\lambda\f$.
///
/// For each grid point the complex variable \f$L = \exp(i\lambda)\f$ is
/// computed and substituted into the series together with the given numeric
/// values of \f$X\f$ and \f$\bar{X}\f$. The series is interpreted as a
/// function of the independent variables \f$X\f$, \f$\bar{X}\f$ and
/// \f$L\f$, with \f$L^{-1}\f$ deriving from \f$L\f$ via \p obake::pow.
///
/// @param s the Keplerian series to evaluate (e.g., the output of
/// \p series::z1 or \p series::z2).
/// @param X_val the numeric value of the variable \f$X\f$.
/// @param Xc_val the numeric value of the variable \f$\bar{X}\f$.
/// @param lambda_start the initial value of the angle \f$\lambda\f$.
/// @param lambda_end the final value of the angle \f$\lambda\f$.
/// @param num_points the number of grid points.
///
/// @return a vector holding the values of \p s at the grid points.
///
/// @throws ::std::invalid_argument if \p num_points is zero.
[[nodiscard]] ::std::vector<cnum_t> evaluate_on_lambda_range(const series &s, cnum_t X_val,
                                                             cnum_t Xc_val,
                                                             real_t lambda_start,
                                                             real_t lambda_end,
                                                             ::std::size_t num_points);

} // namespace kepler

#endif
