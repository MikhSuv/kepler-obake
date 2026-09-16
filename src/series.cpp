// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "series.hpp"

#include <utility>

#include <mp++/rational.hpp>

#include <obake/math/pow.hpp>
#include <obake/math/truncate_p_degree.hpp>
#include <obake/power_series/power_series.hpp>
#include <obake/symbols.hpp>

#include "types.hpp"

namespace kepler
{

namespace detail
{
// TODO:  Написать разложения sin z и cos z

pser_t sqrt_one_minus_t(pser_t t, ::std::int64_t truncation_degree)
{
    pser_t sqrt{rat_t{1}};
    for (::std::int64_t k = 1; k <= truncation_degree; ++k) {
        sqrt += mppp::binomial(rat_t{1, 2}, k) * obake::pow(t, k);
    }
    return sqrt;
}

} // namespace detail

series::series(pser_t s, ::std::int64_t truncation_degree, obake::symbol_set variables)
    : m_series(::std::move(s)), m_truncation_degree(truncation_degree),
      m_variables(::std::move(variables))
{
}

series series::z1(::std::int64_t truncation_degree)
{
    auto [X, Xc, L] = obake::make_p_series<pser_t>("X", "Xc", "L");
    auto factor = I * rat_t{1, 2} * (X * obake::pow(L, -1) - Xc * L);
    auto t = rat_t{-1, 4} * X * Xc;
    auto sqrt = detail::sqrt_one_minus_t(t, truncation_degree);
    auto z1 = factor * sqrt;
    obake::truncate_p_degree(z1, truncation_degree, obake::symbol_set{"X", "Xc"});
    return series{::std::move(z1), truncation_degree, obake::symbol_set{"X", "Xc"}};
}

series series::z2(::std::int64_t truncation_degree)
{
    auto [X, Xc, L] = obake::make_p_series<pser_t>("X", "Xc", "L");
    auto factor = rat_t{1, 2} * (X * obake::pow(L, -1) + Xc * L);
    auto t = rat_t{-1, 4} * X * Xc;
    auto sqrt = detail::sqrt_one_minus_t(t, truncation_degree);
    auto z2 = factor * sqrt;
    obake::truncate_p_degree(z2, truncation_degree, obake::symbol_set{"X", "Xc"});
    return series{::std::move(z2), truncation_degree, obake::symbol_set{"X", "Xc"}};
}

const pser_t &series::get_series() const
{
    return m_series;
}

::std::int64_t series::get_truncation_degree() const
{
    return m_truncation_degree;
}

const obake::symbol_set &series::get_variables() const
{
    return m_variables;
}

} // namespace kepler
