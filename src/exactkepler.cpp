#include <mp++/real.hpp>

#include "kepler/types.hpp"

namespace exact_kepler
{
kepler::real_t z1(kepler::real_t eccentricity, kepler::real_t Mean_anomaly)
{
    return eccentricity * mppp::sin(Mean_anomaly);
}

kepler::real_t z2(kepler::real_t eccentricity, kepler::real_t Mean_anomaly)
{
    return eccentricity * mppp::cos(Mean_anomaly);
}

kepler::real_t z3(kepler::real_t eccentricity, kepler::real_t Mean_anomaly)
{
    kepler::real_t Eccentric_anomaly{Mean_anomaly};

    // Solve E = M + e*sin(E) by fixed-point iteration. The tolerance is tied
    // to the actual precision of the input argument (a few ulps above the
    // 2^-prec rounding noise), rather than to a hard-coded bit count.
    const auto prec = Mean_anomaly.get_prec();
    const auto eps = mppp::exp2(mppp::real{4 - static_cast<long>(prec), prec});

    while (mppp::abs(Eccentric_anomaly - Mean_anomaly -
                     eccentricity * mppp::sin(Eccentric_anomaly)) > eps) {
        Eccentric_anomaly = Mean_anomaly + eccentricity * mppp::sin(Eccentric_anomaly);
    }
    return Eccentric_anomaly - Mean_anomaly;
}

} // namespace exact_kepler
