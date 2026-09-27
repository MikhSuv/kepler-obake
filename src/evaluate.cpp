// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "kepler/evaluate.hpp"

#include <stdexcept>

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
/// \note Touches obake internals (\c _get_s_table() and
/// \c detail::sm_intersect_idx) to walk the coefficient table directly.
cnum_t evaluate_raw(const pser_t &series, const obake::symbol_map<cnum_t> &symbol_map)
{
    const auto &symbol_set = series.get_symbol_set();
    const auto symbol_inersect = obake::detail::sm_intersect_idx(symbol_map, symbol_set);

    cnum_t retval{0, 0};
    for (const auto &tab : series._get_s_table()) {
        for (const auto &t_it : tab) {
            auto k_val = obake::key_evaluate(t_it.first, symbol_inersect, symbol_set);
            const auto &complex_num = t_it.second;
            retval += k_val * cnum_t{complex_num.real(), complex_num.imag()};
        }
    }
    return retval;
}

} // namespace detail

cnum_t evaluate(const series &series, const obake::symbol_map<cnum_t> &symbol_map)
{
    return detail::evaluate_raw(series.get_series(), symbol_map);
}

::std::vector<cnum_t> evaluate_on_lambda_range(const series &series, cnum_t X_val, cnum_t Xc_val,
                                               real_t lambda_start, real_t lambda_end,
                                               ::std::size_t num_points)
{
    if (num_points == 0U) {
        throw ::std::invalid_argument("The number of grid points must be positive");
    }

    ::std::vector<cnum_t> retval;
    retval.reserve(num_points);

    // The distance between consecutive grid points.
    const auto delta = (num_points > 1U)
                           ? (lambda_end - lambda_start) / static_cast<double>(num_points - 1U)
                           : 0.0;

    for (::std::size_t i = 0U; i < num_points; ++i) {
        const auto lambda = lambda_start + static_cast<double>(i) * delta;

        // Compute L = exp(i*lambda).
        cnum_t LMD{0., lambda};
        LMD.exp();

        obake::symbol_map<cnum_t> symbol_map;
        symbol_map.emplace("X", X_val);
        symbol_map.emplace("Xc", Xc_val);
        symbol_map.emplace("L", LMD);

        retval.push_back(evaluate(series, symbol_map));
    }

    return retval;
}

::std::vector<real_t> evaluate_exact_func(real_t (&func)(real_t, real_t),
                                          const real_t &eccentricity, real_t lambda_start,
                                          real_t lambda_end, ::std::size_t num_points)
{
    if (num_points == 0U) {
        throw ::std::invalid_argument("The number of grid points must be positive");
    }
    ::std::vector<real_t> retval;
    const auto delta = (num_points > 1U)
                           ? (lambda_end - lambda_start) / static_cast<double>(num_points - 1U)
                           : 0.0;
    for (::std::size_t i = 0U; i < num_points; i++) {
        const auto lambda = lambda_start + static_cast<double>(i) * delta;
        retval.push_back(func(eccentricity, lambda));
    }

    return retval;
}

} // namespace kepler
