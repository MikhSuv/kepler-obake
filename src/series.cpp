// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "kepler/series.hpp"

#include <fstream>
#include <stdexcept>
#include <utility>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

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

pser_t sqrt_one_minus(const pser_t &base, ::std::int64_t truncation_degree)
{
    // The binomial series of sqrt(1 - t) reads
    //   sqrt(1 - t) = sum_{k>=0} binom(1/2, k) (-t)^k,
    // where the factor (-1)^k is carried by the argument rather than by the
    // coefficient, so that the caller can substitute the already-signed
    // argument. mppp::binomial() returns the coefficient exactly.
    pser_t sqrt_series{rat_t{1}};
    for (::std::int64_t order = 1; order <= truncation_degree; ++order) {
        sqrt_series += ::mppp::binomial(rat_t{1, 2}, order) * ::obake::pow(base, order);
    }
    return sqrt_series;
}

pser_t sin(const pser_t &base, ::std::int64_t truncation_degree)
{
    // Taylor expansion of sin(t) = t - t^3/3! + t^5/5! - ... evaluated with a
    // running term rather than with explicit factorials, so that all the
    // rational arithmetic stays exact.
    pser_t series_term{base};
    pser_t sin_series{series_term};
    for (::std::int64_t order = 1; order <= truncation_degree; ++order) {
        series_term *= rat_t{-1, (2 * order) * ((2 * order) + 1)} * ::obake::pow(base, 2);
        sin_series += series_term;
    }
    return sin_series;
}

pser_t cos(const pser_t &base, ::std::int64_t truncation_degree)
{
    // Taylor expansion of cos(t) = 1 - t^2/2! + t^4/4! - ..., again with a
    // running term.
    pser_t series_term{rat_t{1}};
    pser_t cos_series{series_term};
    for (::std::int64_t order = 1; order <= truncation_degree; ++order) {
        series_term *= rat_t{-1, (2 * order) * ((2 * order) - 1)} * ::obake::pow(base, 2);
        cos_series += series_term;
    }
    return cos_series;
}

} // namespace detail

series::series(pser_t ser, ::std::int64_t truncation_degree, ::obake::symbol_set variables)
    : m_series(::std::move(ser)), m_truncation_degree(truncation_degree),
      m_variables(::std::move(variables))
{
}

series series::z1(::std::int64_t truncation_degree)
{
    // NOTE: NOLINT is required on the structured bindings below. The check
    // applies a hard-coded two-character minimum to binding names and does
    // not consult IgnoredVariableNames for them, so the single-letter
    // Poisson variables cannot be exempted through the configuration.
    // NOLINTNEXTLINE(readability-identifier-length)
    const auto [X, Xc, L] = ::obake::make_p_series<pser_t>("X", "Xc", "L");

    // The harmonic prefactor of z1 carries the factor i/2; z2 below has the
    // same prefactor without it, which is the whole difference between the
    // sine and the cosine combinations of the Poisson variables.
    const pser_t prefactor = I * rat_t{1, 2} * (X * ::obake::pow(L, -1) - Xc * L);

    // The square root is expanded in a dummy symbol t, which is replaced by
    // its value only once the expansion is complete.
    // NOLINTNEXTLINE(readability-identifier-length)
    const auto [t] = ::obake::make_p_series<pser_t>("t");
    ::obake::symbol_map<pser_t> subs_map{{"t", rat_t{-1, 4} * X * Xc}};
    pser_t sqrt_series = ::obake::subs(detail::sqrt_one_minus(t, truncation_degree), subs_map);
    sqrt_series = ::obake::trim(::std::move(sqrt_series));

    // NOTE: the truncation is imposed after the substitution. Truncating the
    // binomial series in t before substituting would cap the degree in the
    // variables at half of what is requested, because t carries degree two in
    // X and Xc.
    ::obake::set_truncation(sqrt_series, truncation_degree, ::obake::symbol_set{"X", "Xc"});

    pser_t z1_ser = prefactor * sqrt_series;
    return series{::std::move(z1_ser), truncation_degree, ::obake::symbol_set{"X", "Xc"}};
}

series series::z2(::std::int64_t truncation_degree)
{
    // NOLINTNEXTLINE(readability-identifier-length)
    const auto [X, Xc, L] = ::obake::make_p_series<pser_t>("X", "Xc", "L");

    // Harmonic prefactor of z2: identical to the one of z1 up to the sign of
    // the imaginary factor.
    const pser_t prefactor = rat_t{1, 2} * (X * ::obake::pow(L, -1) + Xc * L);

    // NOLINTNEXTLINE(readability-identifier-length)
    const auto [t] = ::obake::make_p_series<pser_t>("t");
    ::obake::symbol_map<pser_t> subs_map{{"t", rat_t{-1, 4} * X * Xc}};
    pser_t sqrt_series = ::obake::subs(detail::sqrt_one_minus(t, truncation_degree), subs_map);
    sqrt_series = ::obake::trim(::std::move(sqrt_series));

    // See the note in z1(): the truncation has to follow the substitution.
    ::obake::set_truncation(sqrt_series, truncation_degree, ::obake::symbol_set{"X", "Xc"});

    pser_t z2_ser = prefactor * sqrt_series;
    return series{::std::move(z2_ser), truncation_degree, ::obake::symbol_set{"X", "Xc"}};
}

series series::z3(::std::int64_t truncation_degree)
{
    const auto [Z1, Z2] = ::obake::make_p_series<pser_t>("z1", "z2");
    const ::obake::symbol_set z_vars{"z1", "z2"};

    // Successive approximation of the fixed point of
    //   w = z1 * cos(w) + z2 * sin(w),
    // seeded with the leading-order guess w1 = z1.
    pser_t z3_ser = Z1;

    // NOLINTNEXTLINE(readability-identifier-length)
    const auto [t] = ::obake::make_p_series<pser_t>("t");
    for (::std::int64_t order = 2; order <= truncation_degree; ++order) {

        // Expand the trigonometric functions in the dummy symbol t, which is
        // then replaced by the current approximation. Going through a dummy
        // symbol keeps the Taylor expansion independent of the (possibly
        // very long) running series.
        ::obake::symbol_map<pser_t> subs_map{{"t", z3_ser}};
        pser_t cos_series = ::obake::subs(detail::cos(t, order), subs_map);
        pser_t sin_series = ::obake::subs(detail::sin(t, order), subs_map);

        // NOTE: one order of headroom is granted to the intermediate result.
        // The product with Z1 or Z2 raises the degree by one, so truncating
        // at 'order' here would discard the term that actually carries the
        // requested accuracy.
        ::obake::set_truncation(cos_series, order + 1, z_vars);
        ::obake::set_truncation(sin_series, order + 1, z_vars);

        z3_ser = Z1 * cos_series + Z2 * sin_series;
    }

    z3_ser = ::obake::trim(::std::move(z3_ser));
    ::obake::set_truncation(z3_ser, truncation_degree, z_vars);
    return series{::std::move(z3_ser), truncation_degree, z_vars};
}

series series::z4(::std::int64_t truncation_degree)
{
    const auto [Z3] = ::obake::make_p_series<pser_t>("z3");
    const ::obake::symbol_set z_vars{"z3"};

    pser_t z4_ser = detail::sin(Z3, truncation_degree);
    ::obake::set_truncation(z4_ser, truncation_degree, z_vars);
    return series{::std::move(z4_ser), truncation_degree, z_vars};
}

series series::z5(::std::int64_t truncation_degree)
{
    const auto [Z3] = ::obake::make_p_series<pser_t>("z3");
    const ::obake::symbol_set z_vars{"z3"};

    pser_t z5_ser = detail::cos(Z3, truncation_degree);
    ::obake::set_truncation(z5_ser, truncation_degree, z_vars);
    return series{::std::move(z5_ser), truncation_degree, z_vars};
}

series series::z3_z1z2(const series &z1, const series &z2, const series &z3)
{
    // The truncation degree of the z1/z2 expansions bounds the degree of the
    // final z3, because a monomial of total degree d in (z1, z2) maps to a
    // monomial of total degree d in (X, Xc).
    const ::std::int64_t truncation_degree = z1.get_truncation_degree();
    const ::obake::symbol_set amp_vars{"X", "Xc"};

    // The incoming expansions are re-truncated to the common degree before
    // being substituted, so that no term above the target degree can be
    // produced and later cancelled out.
    pser_t z1_ser = z1.get_series();
    pser_t z2_ser = z2.get_series();
    ::obake::set_truncation(z1_ser, truncation_degree, amp_vars);
    ::obake::set_truncation(z2_ser, truncation_degree, amp_vars);

    // The template is truncated in its own (z1, z2) variables, which is the
    // wrong space once the expansion is closed, hence the truncation policy
    // is dropped before the substitution and reinstalled afterwards.
    pser_t z3_template = z3.get_series();
    ::obake::unset_truncation(z3_template);

    ::obake::symbol_map<pser_t> subs_map{{"z1", z1_ser}, {"z2", z2_ser}};
    pser_t result = ::obake::trim(::obake::subs(z3_template, subs_map));
    ::obake::set_truncation(result, truncation_degree, amp_vars);

    // The symbol set is read before the payload is handed over, so that the
    // move below cannot invalidate it.
    ::obake::symbol_set result_vars = result.get_symbol_set();
    return series{::std::move(result), truncation_degree, ::std::move(result_vars)};
}

const pser_t &series::get_series() const
{
    return m_series;
}

::std::int64_t series::get_truncation_degree() const
{
    return m_truncation_degree;
}

const ::obake::symbol_set &series::get_variables() const
{
    return m_variables;
}

void series::save(const ::std::string &path) const
{
    ::std::ofstream ofs(path);
    if (!ofs) {
        throw ::std::runtime_error("Unable to open file for writing: " + path);
    }
    ::boost::archive::text_oarchive oarchive(ofs);
    oarchive << *this;
    if (!ofs) {
        throw ::std::runtime_error("Error while writing file: " + path);
    }
}

series series::load(const ::std::string &path)
{
    ::std::ifstream ifs(path);
    if (!ifs) {
        throw ::std::runtime_error("Unable to open file for reading: " + path);
    }
    series retval;
    ::boost::archive::text_iarchive iarchive(ifs);
    iarchive >> retval;
    return retval;
}

} // namespace kepler
