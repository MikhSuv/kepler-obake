// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef KEPLER_TYPES_HPP
#define KEPLER_TYPES_HPP

#include <complex>
#include <cstdint>

#include <mp++/complex.hpp>
#include <mp++/rational.hpp>
#include <mp++/real.hpp>

#include <obake/polynomials/packed_monomial.hpp>
#include <obake/power_series/power_series.hpp>

namespace kepler
{

/// Rational number type with fixed 128-bit limb.
///
/// Exact arithmetic is used for the coefficients of the symbolic series: all
/// the operations involved in their construction (binomial coefficients,
/// Taylor recurrences) yield rational numbers, so keeping them exact removes
/// any need for a separate rounding/truncation bookkeeping.
using rat_t = ::mppp::rational<1>;

/// Complex coefficient type with rational real and imaginary parts.
using cf_t = ::std::complex<rat_t>;

/// Arbitrary-precision complex floating-point type used for numerical
/// evaluation of the series.
using cnum_t = ::mppp::complex;

/// Arbitrary-precision real floating-point type.
using real_t = ::mppp::real;

/// Monomial key type: a packed multivariate monomial with 64-bit exponents.
///
/// Negative exponents are allowed and are used to represent inverse powers
/// of the variables, most notably \f$\bar{\Lambda} = \Lambda^{-1}\f$.
using mono_t = ::obake::packed_monomial<::std::int64_t>;

/// Poisson (power) series type with complex rational coefficients.
///
/// The variable order is \f$X, \bar{X}, \Lambda\f$. The variable \f$\Lambda\f$
/// and its powers are treated as parameters by the truncation degree, which
/// is set explicitly through obake's truncation API.
using pser_t = ::obake::p_series<mono_t, cf_t>;

/// The imaginary unit as a complex rational coefficient.
///
/// Defined in terms of rat_t literals, hence the initialisation reduces to
/// mpz_set_si() calls on a zero-sized allocator and cannot allocate memory
/// in practice. NOLINT is needed because the static analyzer cannot prove
/// that mp++'s constructors are non-throwing.
// NOLINTNEXTLINE(bugprone-throwing-static-initialization)
inline const cf_t I{rat_t{0}, rat_t{1}};

} // namespace kepler

#endif
