// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef KEPLER_EXACTKEPLER_HPP
#define KEPLER_EXACTKEPLER_HPP

#include <mp++/real.hpp>

#include "kepler/types.hpp"

namespace kepler::exact_kepler
{

/// Eccentricity times the sine of the mean anomaly, \f$z_1 = e \sin M\f$.
///
/// The arguments are the Keplerian elements in their usual form: the
/// eccentricity @p eccentricity and the mean anomaly @p mean_anomaly
/// (in radians).
///
/// @param eccentricity the orbital eccentricity.
/// @param mean_anomaly the mean anomaly.
///
/// @return the value of \f$e \sin M\f$.
[[nodiscard]] kepler::real_t z1(const kepler::real_t &eccentricity,
                                const kepler::real_t &mean_anomaly);

/// Eccentricity times the cosine of the mean anomaly, \f$z_2 = e \cos M\f$.
///
/// @param eccentricity the orbital eccentricity.
/// @param mean_anomaly the mean anomaly.
///
/// @return the value of \f$e \cos M\f$.
[[nodiscard]] kepler::real_t z2(const kepler::real_t &eccentricity,
                                const kepler::real_t &mean_anomaly);

/// The equation of the centre, \f$z_3 = \nu - M = E - M\f$.
///
/// Equivalently, the difference between the eccentric anomaly @f$E@f$ and
/// the mean anomaly, obtained as the fixed point of Kepler's equation
/// \f[E = M + e \sin E]\f$.
///
/// The equation is solved by plain fixed-point iteration,
/// \f$E \leftarrow M + e \sin E\f$, which converges for every
/// \f$0 \le e < 1\f$. The stopping tolerance is derived from the working
/// precision of @p mean_anomaly rather than from a hard-coded bit count, so
/// that the loop always terminates in a number of steps consistent with the
/// precision the caller asked for.
///
/// @param eccentricity the orbital eccentricity.
/// @param mean_anomaly the mean anomaly.
///
/// @return the value of \f$E - M\f$.
[[nodiscard]] kepler::real_t z3(const kepler::real_t &eccentricity,
                                const kepler::real_t &mean_anomaly);

} // namespace kepler::exact_kepler

#endif
