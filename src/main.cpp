// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

#include <mp++/real.hpp>

#include "evaluate.hpp"
#include "series.hpp"
#include "types.hpp"

int main()
{
    ::std::cout.precision(25);

    auto z1 = kepler::series::z1(1000);
    auto z2 = kepler::series::z2(1000);

    ::std::cout << z1.get_series() << '\n' << '\n';
    ::std::cout << z2.get_series() << '\n' << '\n';

    // Numeric check of z1 at eccentricity e = 0.4 over a full lambda cycle.
    kepler::real_t e{0.4};
    kepler::cnum_t X{mppp::sqrt(kepler::real_t{2}) * mppp::sqrt(1 - mppp::sqrt(1 - e * e)), 0};
    const kepler::real_t PI{mppp::real_pi(200)};

    const auto vals = kepler::evaluate_on_lambda_range(z1, X, X, kepler::real_t{0}, 2 * PI, 10);
    for (const auto &v : vals) {
        ::std::cout << v << '\n';
    }

    return 0;
}
