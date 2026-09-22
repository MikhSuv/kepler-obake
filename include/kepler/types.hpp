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
using rat_t = mppp::rational<1>;

/// Complex coefficient type with rational real and imaginary parts.
using cf_t = ::std::complex<rat_t>;

/// Arbitrary-precision complex floating-point type used for numerical
/// evaluation of the series.
using cnum_t = mppp::complex;

/// Arbitrary-precision real floating-point type.
using real_t = mppp::real;

/// Monomial key type: a packed multivariate monomial with 64-bit exponents.
using mono_t = obake::packed_monomial<::std::int64_t>;

/// Poisson (power) series type with complex rational coefficients.
using pser_t = obake::p_series<mono_t, cf_t>;

/// The imaginary unit as a complex rational coefficient.
inline const cf_t I{rat_t{0}, rat_t{1}};

} // namespace kepler

#endif
