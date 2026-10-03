// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "kepler/exactkepler.hpp"

#include <mp++/real.hpp>

#include "kepler/types.hpp"

namespace kepler::exact_kepler
{

kepler::real_t z1(const kepler::real_t &eccentricity, const kepler::real_t &mean_anomaly)
{
    return eccentricity * ::mppp::sin(mean_anomaly);
}

kepler::real_t z2(const kepler::real_t &eccentricity, const kepler::real_t &mean_anomaly)
{
    return eccentricity * ::mppp::cos(mean_anomaly);
}

kepler::real_t z3(const kepler::real_t &eccentricity, const kepler::real_t &mean_anomaly)
{
    // Kepler's equation is  E = M + e sin(E). The map  E -> M + e sin(E)  is a
    // contraction with modulus at most e, so iterating it converges for every
    // 0 <= e < 1; the fixed point is the eccentric anomaly.
    kepler::real_t ecc_anomaly{mean_anomaly};

    // The stopping tolerance is tied to the actual precision of the input
    // argument (a few ulps above the 2^-prec rounding noise), rather than to
    // a hard-coded bit count, so that the loop always terminates in a number
    // of steps consistent with the precision the caller asked for.
    const auto prec = mean_anomaly.get_prec();
    const auto eps = ::mppp::exp2(::mppp::real{4 - static_cast<long>(prec), prec});

    // The residual is measured on the left-hand side of Kepler's equation, so
    // the loop stops as soon as the two sides agree to within eps.
    while (::mppp::abs(ecc_anomaly - mean_anomaly - eccentricity * ::mppp::sin(ecc_anomaly))
           > eps) {
        ecc_anomaly = mean_anomaly + eccentricity * ::mppp::sin(ecc_anomaly);
    }

    // Return the equation of the centre rather than the eccentric anomaly.
    return ecc_anomaly - mean_anomaly;
}

kepler::real_t z4(const kepler::real_t &eccentricity, const kepler::real_t &mean_anomaly)
{
    return ::mppp::sin(z3(eccentricity, mean_anomaly));
}

kepler::real_t z5(const kepler::real_t &eccentricity, const kepler::real_t &mean_anomaly)
{
    return ::mppp::cos(z3(eccentricity, mean_anomaly));
}

} // namespace kepler::exact_kepler
