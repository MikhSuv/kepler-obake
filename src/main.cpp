#include <complex>
#include <cstdint>
#include <iostream>
#include <mp++/complex.hpp>
#include <mp++/rational.hpp>
#include <obake/key/key_evaluate.hpp>
#include <obake/math/pow.hpp>
#include <obake/math/truncate_p_degree.hpp>
#include <obake/polynomials/packed_monomial.hpp>
#include <obake/power_series/power_series.hpp>
#include <obake/symbols.hpp>

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
  auto z1 = z1_series(50);
  auto z2 = z2_series(10);
  std::cout << z1 << '\n' << '\n';
  std::cout << z2 << '\n';

  return 0;
}
