#include <complex>
#include <iostream>
#include <mp++/complex.hpp>
#include <mp++/rational.hpp>
#include <obake/key/key_evaluate.hpp>
#include <obake/math/pow.hpp>
#include <obake/math/truncate_p_degree.hpp>
#include <obake/polynomials/packed_monomial.hpp>
#include <obake/polynomials/polynomial.hpp>
#include <obake/power_series/power_series.hpp>
#include <obake/symbols.hpp>

using rat_t = mppp::rational<1>;
using cf_t = std::complex<rat_t>;
using cnum_t = mppp::complex;
using mono_t = obake::packed_monomial<std::int64_t>;
using poly_t = obake::p_series<mono_t, cf_t>;

template <typename Series, typename Map>
cnum_t my_evaluate(const Series &s, const Map &sm) {
  const auto &ss = s.get_symbol_set();
  const auto si = obake::detail::sm_intersect_idx(sm, ss);

  cnum_t retval{0, 0};
  for (const auto &t : s) {
    auto k_val = obake::key_evaluate(t.first, si, ss);
    const auto &c = t.second;
    retval += k_val * cnum_t{c.real(), c.imag()};
  }
  return retval;
}

int main() {
  cf_t i{rat_t{0}, rat_t{1}};

  auto [a, b, c] = obake::make_p_series<poly_t>("a", "b", "c");

  auto base = a * obake::pow(a, -1) + obake::pow(b, 2) + obake::pow(c, -1);
  auto result = obake::pow(base, 5u);

  obake::symbol_map<cnum_t> num_map{
      {"a", cnum_t{5, 0}}, {"b", cnum_t{0.5, 0}}, {"c", cnum_t{1.333, 0}}};

  cnum_t num = my_evaluate(result, num_map);

  std::cout << num << '\n';

  return 0;
}
