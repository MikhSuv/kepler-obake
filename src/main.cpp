#include <complex>
#include <cstdint>
#include <iostream>
#include <mp++/complex.hpp>
#include <mp++/rational.hpp>
#include <mp++/real.hpp>
#include <obake/key/key_evaluate.hpp>
#include <obake/math/pow.hpp>
#include <obake/math/truncate_p_degree.hpp>
#include <obake/polynomials/packed_monomial.hpp>
#include <obake/power_series/power_series.hpp>
#include <obake/symbols.hpp>
#include <stdexcept>
#include <vector>

using rat_t = mppp::rational<1>;
using cf_t = std::complex<rat_t>;
using cnum_t = mppp::complex;
using mono_t = obake::packed_monomial<std::int64_t>;
using poly_t = obake::p_series<mono_t, cf_t>;

// Обертка для evaluate, чтобы преобразовать к комплексным числам с
// плавающей точкой.
template <typename Series>
cnum_t my_evaluate(const Series &s, const obake::symbol_map<cnum_t> &sm) {
  const auto &ss = s.get_symbol_set();
  const auto si = obake::detail::sm_intersect_idx(sm, ss);

  cnum_t retval{0, 0};
  for (const auto &tab : s._get_s_table()) {
    for (const auto &t : tab) {
      auto k_val = obake::key_evaluate(t.first, si, ss);
      const auto &c = t.second;
      retval += k_val * cnum_t{c.real(), c.imag()};
    }
  }
  return retval;
}

// Мнимая единица
const cf_t I{rat_t{0}, rat_t{1}};
const mppp::real PI{mppp::real_pi(200)};

/**
 * @brief Evaluate a Poisson series on a uniform grid of the angle
 * \f$\lambda\f$.
 *
 * For each grid point the complex variable \f$L = \exp(i\lambda)\f$ is
 * computed and substituted into the series together with the given numeric
 * values of \f$X\f$ and \f$\bar{X}\f$. The series is interpreted as a
 * function of the independent variables \f$X\f$, \f$\bar{X}\f$ and
 * \f$L\f$, with \f$L^{-1}\f$ deriving from \f$L\f$ via \p obake::pow.
 *
 * @param series the Poisson series to evaluate (e.g., the output of
 * \p z1_series or \p z2_series).
 * @param X_val the numeric value of the variable \f$X\f$.
 * @param Xc_val the numeric value of the variable \f$\bar{X}\f$.
 * @param lambda_start the initial value of the angle \f$\lambda\f$.
 * @param lambda_end the final value of the angle \f$\lambda\f$.
 * @param num_points the number of grid points.
 *
 * @return a vector holding the values of \p series at the grid points.
 *
 * @throws ::std::invalid_argument if \p num_points is zero.
 */
::std::vector<cnum_t> evaluate_on_lambda_range(const poly_t &series,
                                               cnum_t X_val, cnum_t Xc_val,
                                               mppp::real lambda_start,
                                               mppp::real lambda_end,
                                               ::std::size_t num_points) {
  if (num_points == 0) {
    throw ::std::invalid_argument("The number of grid points must be positive");
  }

  ::std::vector<cnum_t> retval;
  retval.reserve(num_points);

  // The distance between consecutive grid points.
  const auto delta =
      (num_points > 1u)
          ? (lambda_end - lambda_start) / static_cast<double>(num_points - 1u)
          : 0.0;

  for (::std::size_t i = 0; i < num_points; ++i) {
    const auto lambda = lambda_start + static_cast<double>(i) * delta;

    // Compute L = exp(i*lambda).
    cnum_t L{0., lambda};
    L.exp();

    obake::symbol_map<cnum_t> sm;
    sm.emplace("X", X_val);
    sm.emplace("Xc", Xc_val);
    sm.emplace("L", L);

    // retval.push_back(my_evaluate(series, sm));
    retval.push_back(my_evaluate(series, sm));
  }

  return retval;
}

poly_t z1_series(int64_t total_degree) {
  auto [X, Xc, L] = obake::make_p_series<poly_t>("X", "Xc", "L");
  auto factor = I * rat_t{0.5} * (X * obake::pow(L, -1) - Xc * L);
  auto tmp = rat_t{-1, 4} * X * Xc;
  poly_t sqrt{rat_t{1}};
  for (int k = 1; k <= total_degree; ++k) {
    sqrt += mppp::binomial(rat_t{1, 2}, k) * obake::pow(tmp, k);
  }
  auto z1 = factor * sqrt;
  obake::truncate_p_degree(z1, total_degree, obake::symbol_set{"X", "Xc"});
  return z1;
}

poly_t z2_series(int64_t total_degree) {
  auto [X, Xc, L] = obake::make_p_series<poly_t>("X", "Xc", "L");
  auto factor = rat_t{0.5} * (X * obake::pow(L, -1) + Xc * L);
  auto tmp = rat_t{-1, 4} * X * Xc;
  poly_t sqrt{rat_t{1}};
  for (int k = 1; k <= total_degree; ++k) {
    sqrt += mppp::binomial(rat_t{1, 2}, k) * obake::pow(tmp, k);
  }
  auto z2 = factor * sqrt;
  obake::truncate_p_degree(z2, total_degree, obake::symbol_set{"X", "Xc"});
  return z2;
}

int main(void) {
  std::cout.precision(25);
  auto z1 = z1_series(1000);
  auto z2 = z2_series(1000);
  std::cout << z1 << '\n' << '\n';
  std::cout << z2 << '\n' << '\n';
  // std::cout << z2 << '\n';
  mppp::real e{0.4};
  cnum_t X{mppp::sqrt(mppp::real{2}) * mppp::sqrt(1 - mppp::sqrt(1 - e * e)),
           0};

  // Evaluate z1 on a grid of lambda angles with fixed X and Xc.
  const auto vals =
      evaluate_on_lambda_range(z1, X, X, mppp::real{0}, 2 * PI, 10);
  for (const auto &v : vals) {
    // std::cout << v.get_real_imag().first << '\n';
    std::cout << v << '\n';
  }

  return 0;
}
