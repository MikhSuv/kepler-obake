#include "kepler/series.hpp"
#include "kepler/types.hpp"
#include <cstdint>
#include <iostream>

int main()
{
    auto z1 = kepler::series::load("z1.txt");
    auto z2 = kepler::series::load("z2.txt");
    std::cout << "Loaded z1 and z2 \n";

    int64_t degree{std::min(z1.get_truncation_degree(), z2.get_truncation_degree())};

    auto z3_in_z1_z2 = kepler::series::z3(degree);
    z3_in_z1_z2.save("z3_in_z1_z2.txt");
    std::cout << "Saved z3 in terms of z1 and z2. Term number:" << z3_in_z1_z2.get_series().size()
              << '\n';

    auto z3 = kepler::series::z3_z1z2(z1, z2, z3_in_z1_z2);
    z3.save("z3.txt");
    std::cout << "Saved z3 in terms of X and Xc. Term number:" << z3.get_series().size() << '\n';
    return 0;
}
