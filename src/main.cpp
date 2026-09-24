// Project: kepler-obake
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "kepler/evaluate.hpp"
#include "kepler/exactkepler.hpp"
#include "kepler/series.hpp"
#include "kepler/types.hpp"
#include <algorithm>
#include <iostream>
#include <mp++/real.hpp>
#include <obake/math/truncate_p_degree.hpp>

int main()
{
    ::std::cout.precision(16);
    //
    int deg = 10;
    auto z1 = kepler::series::z1(deg);
    auto z2 = kepler::series::z2(deg);
    //
    ::std::cout << z1.get_series().size() << '\n' << '\n';
    ::std::cout << z2.get_series().size() << '\n' << '\n';

    // Working precision (in bits) for all the numeric inputs below. mp++
    // deduces the precision of a real from the type of its source (32 bits
    // for int, 53 bits for double), so explicit precision is required to
    // avoid capping the achievable accuracy regardless of the series degree.
    const long wp = 256;
    kepler::real_t e{"0.3", wp};
    kepler::cnum_t X{
        mppp::sqrt(kepler::real_t{2, wp}) *
            mppp::sqrt(kepler::real_t{1, wp} - mppp::sqrt(kepler::real_t{1, wp} - e * e)),
        kepler::real_t{0, wp}};
    const kepler::real_t PI{mppp::real_pi(wp)};

    // const auto vals = kepler::evaluate_on_lambda_range(z1, X, X, kepler::real_t{0}, 2 * PI, 10);
    // for (const auto &v : vals) {
    //     ::std::cout << v << '\n';
    // }

    auto z3 = kepler::series::z3(deg);
    ::std::cout << "Создан ряд z3. количество слагаемых ";
    ::std::cout << z3.get_series().size() << '\n';

    z3 = kepler::series::z3_z1z2(z1, z2, z3);
    ::std::cout << "Создан ряд z3.";
    ::std::cout << z3.get_series() << '\n';

    ::std::size_t N = 1000u;

    const auto vals = kepler::evaluate_on_lambda_range(z3, X, X, kepler::real_t{0, wp}, 2 * PI, N);
    ::std::cout << "Все значения для ряда посчитаны \n";

    const auto exact_vals =
        kepler::evaluate_exact_func(exact_kepler::z3, e, kepler::real_t{0, wp}, 2 * PI, N);
    ::std::cout << "Все значения по формуле посчитаны \n";
    ::std::vector<kepler::real_t> abs_err;
    for (::std::size_t i = 0u; i < N; ++i) {
        abs_err.push_back(mppp::abs(vals[i].get_real_imag().first - exact_vals[i]));
        ::std::cout << vals[i].get_real_imag().first << ' ' << exact_vals[i] << ' ' << abs_err[i]
                    << '\n';
    }
    ::std::cout << "Max error :" << *std::max_element(abs_err.begin(), abs_err.end()) << '\n';
    return 0;
}
