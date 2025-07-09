#ifndef TEST_H
#define TEST_H

#include <cancelator/fuelpin_approximate_cancelator.hpp>

void test(){
    std::cout << "---------------------------------" << std::endl;
    Position low(-75.2626, -75.2626, -110.4155);
    Position high(75.2626, 75.2626, 89.5845);

    std::vector<double> assembly_pitch = {21.5036, 21.5036, 200.};
    std::vector<std::uint32_t> assembly_shape = {7, 7, 1};
    std::vector<double> fuelpin_pitch = {1.2598, 1.2598, 200.};
    std::vector<std::uint32_t> fuelpins_per_assembly = {17, 17, 1};
    std::vector<double> inter_assembly_gap = {0.087, 0.087, 0.};
    const double radius = 0.4058;
    const std::size_t Nr_ = 4;

    FuelPinApproxCancelator fp_cancel(low, high, 
                                      assembly_pitch, assembly_shape,
                                      fuelpin_pitch, fuelpins_per_assembly,
                                      inter_assembly_gap, radius,
                                      Nr_, FuelPinApproxCancelator::Orientation::Z);

    BankedParticle p;
    p.r = Position(51.8255, 34.145,-10.);
    bool t = fp_cancel.add_particle(p);
    std::cout << " <<< " << t << std::endl;

    std::cout << "---------------------------------" << std::endl;
}

#endif