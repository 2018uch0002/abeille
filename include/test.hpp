#ifndef TEST_H
#define TEST_H

#include <tallies/cylinder_filter.hpp>

void test() {
  std::cout << "We are inside of the test function. \n\n" << std::endl;

  Position origin(-10.08, -10.08, -10.71);
  const double radius = 0.54;
  const double pitch = 1.26;
  const double dz = 2.142;
  CylinderFilter cy_filter(origin, radius, pitch, pitch, dz, 17, 17, 10,
                           CylinderFilter::Orientation::Z, 1);

  double d_flight =  10.119985517406355;
  Position r(1.08006565031154, 1.26,-3.16518597623835);
  Direction u(-0.555535221291904,-0.534165667308047,0.637218689127606);
  
  Tracker trkr(r, u, true);
  auto index_distance = cy_filter.get_indices_tracklength_with_position(trkr, d_flight);
  
  std::cout << "length of the vector = " << index_distance.size() << std::endl;
  for (auto& p : index_distance) {
    std::cout << "index: " << p.index[0] << ", " << p.index[1] << ", "
              << p.index[2] << "\tdistance = " << p.distance << "\tr0 = " << p.r0 << std::endl;
  }


  std::cout << "\n\nWe are exiting the test function. \n\n" << std::endl;
}

#endif