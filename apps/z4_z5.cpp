#include "kepler/series.hpp"
#include <iostream>

int main()
{
    auto z3 = kepler::series::load("z3.txt");
    auto degree = z3.get_truncation_degree();
    auto z4_in_z3 = kepler::series::z4(degree);
    z4_in_z3.save("z4_in_z3.txt");
    auto z5_in_z3 = kepler::series::z5(degree);
    z5_in_z3.save("z5_in_z3.txt");

    auto z4 = kepler::series::z4z5_z3(z3, z4_in_z3);
    z4.save("z4.txt");
    std::cout << "z4: Number of terms: " << z4.get_series().size() << '\n';
    auto z5 = kepler::series::z4z5_z3(z3, z5_in_z3);
    z5.save("z5.txt");
    std::cout << "z5: Number of terms: " << z5.get_series().size() << '\n';

    return 0;
}
