#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <iostream>

namespace bmp = boost::multiprecision;

using rat_t   = bmp::mpq_rational;                                        // rational
using real_t  = bmp::mpfr_float_100;                                      // arbitrary precision float
using cf_t    = std::complex<rat_t>;                                      // complex rational
using cnum_t  = bmp::number<bmp::backends::mpc_complex_backend<50>>;      // complex float (MPC)

int main()
{
    // --- Exact rational complex arithmetic ---
    cf_t z1{rat_t{1, 4}, rat_t{2, 5}};
    cf_t z2{rat_t{3, 7}, rat_t{-1, 3}};
    auto prod = z1 * z2;
    std::cout << "z1 * z2 = (" << prod.real() << ", " << prod.imag() << ")i\n";

    // --- Arbitrary precision complex float arithmetic ---
    cnum_t x(real_t{"1.12345678901234567890123456789"},
             real_t{"2.98765432109876543210987654321"});
    cnum_t y(real_t{"0.55555555555555555555555555555"},
             real_t{"-1.11111111111111111111111111111"});
    std::cout << "x + y = " << x + y << '\n';
    std::cout << "x * y = " << x * y << '\n';

    // --- Converting rational complex -> float complex ---
    cnum_t fz(real_t{z1.real()}, real_t{z1.imag()});
    std::cout << "z1 as float = " << fz << '\n';
}