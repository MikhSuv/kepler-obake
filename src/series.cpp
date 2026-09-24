// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "kepler/series.hpp"

#include <utility>

#include <mp++/rational.hpp>

#include <obake/math/pow.hpp>
#include <obake/math/subs.hpp>
#include <obake/math/trim.hpp>
#include <obake/math/truncate_p_degree.hpp>
#include <obake/power_series/power_series.hpp>
#include <obake/symbols.hpp>

#include "kepler/types.hpp"

namespace kepler
{

namespace detail
{

pser_t sqrt_one_minus_t(const pser_t &t, ::std::int64_t truncation_degree)
{
    pser_t sqrt{rat_t{1}};
    for (::std::int64_t k = 1; k <= truncation_degree; ++k) {
        sqrt += mppp::binomial(rat_t{1, 2}, k) * obake::pow(t, k);
    }
    return sqrt;
}

pser_t sin(const pser_t &t, ::std::int64_t truncation_degree)
{
    pser_t series_term{t};
    pser_t sin_series{series_term};
    for (::std::int64_t k = 1; k <= truncation_degree; ++k) {
        series_term *= rat_t{-1, (2 * k) * (2 * k + 1)} * obake::pow(t, 2);
        sin_series += series_term;
    }
    return sin_series;
}

pser_t cos(const pser_t &t, ::std::int64_t truncation_degree)
{
    pser_t series_term{rat_t{1}};
    pser_t cos_series = series_term;
    for (::std::int64_t k = 1; k <= truncation_degree; ++k) {
        series_term *= rat_t{-1, (2 * k) * (2 * k - 1)} * obake::pow(t, 2);
        cos_series += series_term;
    }
    return cos_series;
}

} // namespace detail

series::series(const pser_t &s, ::std::int64_t truncation_degree,
               const obake::symbol_set &variables)
    : m_series(::std::move(s)), m_truncation_degree(truncation_degree),
      m_variables(::std::move(variables))
{
}
series series::z1(::std::int64_t truncation_degree)
{
    auto [X, Xc, L] = obake::make_p_series<pser_t>("X", "Xc", "L");
    auto factor = I * rat_t{1, 2} * (X * obake::pow(L, -1) - Xc * L);
    auto t = rat_t{-1, 4} * X * Xc;
    // TODO: Переделать через подстановку
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
    // TODO: Переделать через подстановку
    auto sqrt = detail::sqrt_one_minus_t(t, truncation_degree);
    auto z2 = factor * sqrt;
    obake::truncate_p_degree(z2, truncation_degree, obake::symbol_set{"X", "Xc"});
    return series{::std::move(z2), truncation_degree, obake::symbol_set{"X", "Xc"}};
}

series series::z3(::std::int64_t truncation_degree)
{
    auto [Z1, Z2] = obake::make_p_series<pser_t>("z1", "z2");
    pser_t w = Z1;
    auto [t] = obake::make_p_series<pser_t>("t");
    obake::symbol_set ss{"z1", "z2"};
    pser_t cos{};
    pser_t sin{};

    for (::std::int64_t n = 2; n <= truncation_degree; ++n) {

        obake::symbol_map<pser_t> sm{{"t", w}};
        cos = detail::cos(t, n);
        sin = detail::sin(t, n);
        cos = obake::subs(cos, sm);
        sin = obake::subs(sin, sm);
        w = Z1 * cos + Z2 * sin;

        obake::truncate_p_degree(w, n + 1, ss);
    }

    obake::truncate_p_degree(w, truncation_degree, ss);

    return series{std::move(w), truncation_degree, ss};
}

series series::z3_z1z2(const series &z1, const series &z2, const series &z3)
{
    ::std::int64_t truncation_degree = z1.get_truncation_degree();
    obake::symbol_map<pser_t> sm{{"z1", z1.get_series()}, {"z2", z2.get_series()}};
    pser_t res = obake::subs(z3.get_series(), sm);
    res = obake::trim(res);
    obake::symbol_set z3_ss = res.get_symbol_set();
    obake::truncate_p_degree(res, truncation_degree, {"X", "Xc"});
    return series{std::move(res), truncation_degree, z3_ss};
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
