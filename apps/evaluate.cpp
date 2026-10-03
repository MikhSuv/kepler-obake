#include "kepler/evaluate.hpp"
#include "kepler/exactkepler.hpp"
#include "kepler/series.hpp"
#include "kepler/types.hpp"
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

void tabulate(const kepler::series &z_series, auto &&exact_func, const std::string &filename);

int main()
{
    ::std::cout.precision(16);
    auto z1 = kepler::series::load("z1.txt");
    auto z2 = kepler::series::load("z2.txt");
    auto z3 = kepler::series::load("z3.txt");
    auto z4 = kepler::series::load("z4.txt");
    auto z5 = kepler::series::load("z5.txt");

    tabulate(z1, kepler::exact_kepler::z1, "tabulate_z1.txt");
    tabulate(z2, kepler::exact_kepler::z2, "tabulate_z2.txt");
    tabulate(z3, kepler::exact_kepler::z3, "tabulate_z3.txt");
    tabulate(z4, kepler::exact_kepler::z4, "tabulate_z4.txt");
    tabulate(z5, kepler::exact_kepler::z5, "tabulate_z5.txt");
    return 0;
}

void tabulate(const kepler::series &z_series, auto &&exact_func, const std::string &filename)
{
    const long work_precision = 400;
    const size_t num_points = 1000;
    kepler::real_t e{"0.3", work_precision};
    kepler::cnum_t X{mppp::sqrt(kepler::real_t{2, work_precision})
                         * mppp::sqrt(kepler::real_t{1, work_precision}
                                      - mppp::sqrt(kepler::real_t{1, work_precision} - e * e)),
                     kepler::real_t{0, work_precision}};
    const kepler::real_t PI{mppp::real_pi(work_precision)};

    const auto vals = kepler::evaluate_on_lambda_range(
        z_series, X, X, kepler::real_t{0, work_precision}, 2 * PI, num_points);
    const auto exact_vals = kepler::evaluate_exact_func(
        exact_func, e, kepler::real_t{0, work_precision}, 2 * PI, num_points);

    std::ofstream out(filename);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << "# series \t exact_value \t error\n";
    std::vector<kepler::real_t> abs_err;
    for (size_t it = 0; it < num_points; ++it) {
        abs_err.push_back(mppp::abs(vals.at(it).get_real_imag().first - exact_vals.at(it)));
        out << vals.at(it).get_real_imag().first << exact_vals.at(it) << abs_err.at(it) << '\n';
    }
    std::cout << filename << " Max error: " << *std::ranges::max_element(abs_err) << '\n';
    out.close();
}
