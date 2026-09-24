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

    // Solve E = M + e*sin(E) by fixed-point iteration. The tolerance is tied
    // to the actual precision of the input argument (a few ulps above the
    // 2^-prec rounding noise), rather than to a hard-coded bit count.
    const auto prec = M.get_prec();
    const auto eps = mppp::exp2(mppp::real{4 - static_cast<long>(prec), prec});

    do {
        E = M + e * mppp::sin(E);

    } while (mppp::abs(E - M - e * mppp::sin(E)) > eps);
    return E - M;
}

} // namespace exact_kepler
