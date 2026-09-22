// #include <mp++/complex.hpp>
#include <mp++/real.hpp>

#include "kepler/types.hpp"

namespace exact_kepler
{
kepler::real_t z1(kepler::real_t e, kepler::real_t M)
{
    return e * mppp::sin(M);
}

kepler::real_t z2(kepler::real_t e, kepler::real_t M)
{
    return e * mppp::cos(M);
}

kepler::real_t z3(kepler::real_t e, kepler::real_t M)
{
    kepler::real_t E{M};
    const auto eps = mppp::exp2(kepler::real_t{1} - static_cast<kepler::real_t>(128));
    do {
        E = M + e * mppp::sin(E);

    } while (mppp::abs(E - M - e * mppp::sin(E)) > eps);
    return E - M;
}

} // namespace exact_kepler
