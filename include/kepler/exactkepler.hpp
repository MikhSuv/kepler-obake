#ifndef EXACTKEPLER_HPP
#define EXACTKEPLER_HPP
#include <mp++/real.hpp>

#include "types.hpp"

namespace exact_kepler
{
kepler::real_t z1(kepler::real_t e, kepler::real_t M);

kepler::real_t z2(kepler::real_t e, kepler::real_t M);

kepler::real_t z3(kepler::real_t e, kepler::real_t M);

} // namespace exact_kepler
//
#endif
