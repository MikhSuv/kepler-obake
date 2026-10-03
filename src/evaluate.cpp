// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "kepler/evaluate.hpp"

#include <stdexcept>
#include <vector>

#include <mp++/complex.hpp>
#include <mp++/real.hpp>

#include <obake/key/key_evaluate.hpp>
#include <obake/symbols.hpp>

#include "kepler/series.hpp"
#include "kepler/types.hpp"

namespace kepler
{

namespace detail
{

/// Evaluate the underlying obake series, converting complex rational
/// coefficients to the arbitrary-precision complex type \ref cnum_t.
///
/// Rather than going through obake's generic stream inserter, the
/// coefficient table is walked directly and each term is evaluated
/// individually. This keeps the conversion of the rational coefficients to
/// \ref cnum_t exact and lets the arbitrary-precision accumulation happen
/// only once per term.
///
/// @note Touches obake internals (\c _get_s_table() and
/// \c detail::sm_intersect_idx) to walk the coefficient table directly.
cnum_t evaluate_raw(const pser_t &ser, const ::obake::symbol_map<cnum_t> &symbol_map)
{
    const auto &symbol_set = ser.get_symbol_set();
    const auto symbol_intersect = ::obake::detail::sm_intersect_idx(symbol_map, symbol_set);

    cnum_t retval{0, 0};
    for (const auto &tab : ser._get_s_table()) {
        for (const auto &term : tab) {
            const auto key_val = ::obake::key_evaluate(term.first, symbol_intersect, symbol_set);
            const auto &coefficient = term.second;
            retval += key_val * cnum_t{coefficient.real(), coefficient.imag()};
        }
    }
    return retval;
}

} // namespace detail

cnum_t evaluate(const series &ser, const ::obake::symbol_map<cnum_t> &symbol_map)
{
    return detail::evaluate_raw(ser.get_series(), symbol_map);
}

::std::vector<cnum_t> evaluate_on_lambda_range(const series &ser, cnum_t X_val, cnum_t Xc_val,
                                               const real_t &lambda_start, const real_t &lambda_end,
                                               ::std::size_t num_points)
{
    if (num_points == 0U) {
        throw ::std::invalid_argument("The number of grid points must be positive");
    }

    ::std::vector<cnum_t> retval;
    retval.reserve(num_points);

    // The distance between consecutive grid points.
    // NOTE: the divisor is narrowed to double, which caps the resolution of
    // the grid spacing at the precision of a binary64 number. This is
    // intentional: lambda itself is a real_t, so the accuracy of the grid is
    // limited by this conversion rather than by the working precision of the
    // surrounding arithmetic.
    const auto delta = (num_points > 1U)
                           ? (lambda_end - lambda_start) / static_cast<double>(num_points - 1U)
                           : 0.0;

    for (::std::size_t idx = 0U; idx < num_points; ++idx) {
        const auto lambda = lambda_start + static_cast<double>(idx) * delta;

        // Compute L = exp(i*lambda), the Poisson variable standing for the
        // mean longitude.
        cnum_t lmd{0., lambda};
        lmd.exp();

        ::obake::symbol_map<cnum_t> sym_map;
        sym_map.emplace("X", X_val);
        sym_map.emplace("Xc", Xc_val);
        sym_map.emplace("L", lmd);

        retval.push_back(evaluate(ser, sym_map));
    }

    return retval;
}

::std::vector<real_t> evaluate_exact_func(real_t (&func)(const real_t &, const real_t &),
                                          const real_t &eccentricity, const real_t &lambda_start,
                                          const real_t &lambda_end, ::std::size_t num_points)
{
    if (num_points == 0U) {
        throw ::std::invalid_argument("The number of grid points must be positive");
    }

    ::std::vector<real_t> retval;
    retval.reserve(num_points);

    // Same grid, and hence the same double-narrowed spacing, as in
    // evaluate_on_lambda_range().
    const auto delta = (num_points > 1U)
                           ? (lambda_end - lambda_start) / static_cast<double>(num_points - 1U)
                           : 0.0;

    for (::std::size_t idx = 0U; idx < num_points; ++idx) {
        const auto lambda = lambda_start + static_cast<double>(idx) * delta;
        retval.push_back(func(eccentricity, lambda));
    }

    return retval;
}

} // namespace kepler
