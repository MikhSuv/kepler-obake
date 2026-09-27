#ifndef EXACTKEPLER_HPP
#define EXACTKEPLER_HPP
#include <mp++/real.hpp>

#include "types.hpp"

namespace exact_kepler
{
kepler::real_t z1(kepler::real_t eccentricity, kepler::real_t Mean_anomaly);

kepler::real_t z2(kepler::real_t eccentricity, kepler::real_t Mean_anomaly);

kepler::real_t z3(kepler::real_t eccentricity, kepler::real_t Mean_anomaly);

} // namespace exact_kepler
//
#endif
