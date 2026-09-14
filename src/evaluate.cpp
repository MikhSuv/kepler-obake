// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "evaluate.hpp"

#include <stdexcept>

#include <mp++/complex.hpp>
#include <mp++/real.hpp>

#include <obake/key/key_evaluate.hpp>
#include <obake/symbols.hpp>

#include "series.hpp"
#include "types.hpp"

namespace kepler
{

namespace detail
{

/// Evaluate the underlying obake series, converting complex rational
/// coefficients to the arbitrary-precision complex type \ref cnum_t.
///
/// \note Touches obake internals (\c _get_s_table() and
/// \c detail::sm_intersect_idx) to walk the coefficient table directly.
cnum_t evaluate_raw(const pser_t &s, const obake::symbol_map<cnum_t> &sm)
{
    const auto &ss = s.get_symbol_set();
    const auto si = obake::detail::sm_intersect_idx(sm, ss);

    cnum_t retval{0, 0};
    for (const auto &tab : s._get_s_table()) {
        for (const auto &t : tab) {
            auto k_val = obake::key_evaluate(t.first, si, ss);
            const auto &c = t.second;
            retval += k_val * cnum_t{c.real(), c.imag()};
        }
    }
    return retval;
}

} // namespace detail

cnum_t evaluate(const series &s, const obake::symbol_map<cnum_t> &sm)
{
    return detail::evaluate_raw(s.get_series(), sm);
}

::std::vector<cnum_t> evaluate_on_lambda_range(const series &s, cnum_t X_val, cnum_t Xc_val,
                                               real_t lambda_start, real_t lambda_end,
                                               ::std::size_t num_points)
{
    if (num_points == 0u) {
        throw ::std::invalid_argument("The number of grid points must be positive");
    }

    ::std::vector<cnum_t> retval;
    retval.reserve(num_points);

    // The distance between consecutive grid points.
    const auto delta = (num_points > 1u)
                           ? (lambda_end - lambda_start) / static_cast<double>(num_points - 1u)
                           : 0.0;

    for (::std::size_t i = 0u; i < num_points; ++i) {
        const auto lambda = lambda_start + static_cast<double>(i) * delta;

        // Compute L = exp(i*lambda).
        cnum_t L{0., lambda};
        L.exp();

        obake::symbol_map<cnum_t> sm;
        sm.emplace("X", X_val);
        sm.emplace("Xc", Xc_val);
        sm.emplace("L", L);

        retval.push_back(evaluate(s, sm));
    }

    return retval;
}

} // namespace kepler
