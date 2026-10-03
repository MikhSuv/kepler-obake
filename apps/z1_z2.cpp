#include "kepler/series.hpp"
#include "kepler/types.hpp"
#include <iostream>

int main()
{
    constexpr int degree = 100;
    auto z1 = kepler::series::z1(degree);
    std::cout << "z1 series has been created, term number: " << z1.get_series().size() << '\n';
    auto z2 = kepler::series::z2(degree);
    std::cout << "z2 series has been created, term number: " << z2.get_series().size() << '\n';

    z1.save("z1.txt");
    z2.save("z2.txt");
    std::cout << "z1 and z2 are saved";

    return 0;
}
