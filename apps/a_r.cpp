#include <iostream>
#include <kepler/series.hpp>

int main()
{
    auto z1 = kepler::series::load("z1.txt");
    auto z2 = kepler::series::load("z2.txt");
    auto z4 = kepler::series::load("z4.txt");
    auto z5 = kepler::series::load("z5.txt");

    const auto degree = z1.get_truncation_degree();

    {
        auto ra = kepler::series::r_over_a(degree);
        auto rax = kepler::series::r_a_x(z1, z2, z4, z5, ra);
        std::cout << "r/a " << rax.get_series().size() << '\n';
        std::cout << rax.get_series() << '\n';
        rax.save("r_a.txt");
    }
    {
        auto ar = kepler::series::a_over_r(degree);
        auto arx = kepler::series::r_a_x(z1, z2, z4, z5, ar);
        std::cout << "a/r " << arx.get_series().size() << '\n';
        std::cout << arx.get_series() << '\n';
        arx.save("a_r.txt");
    }
    return 0;
}
