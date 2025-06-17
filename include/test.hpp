#ifndef TEST_H
#define TEST_H

#include <tallies/cylinder_filter.hpp>

void test() {
  std::cout << "We are inside of the test function. \n\n" << std::endl;
  Position origin(-10.0784, -10.0784, -100.);
  const double radius = 0.4750;
  const double pitch = 1.2598;
  const double dz = 200.;

  CylinderFilter cy_filter(origin, radius, pitch, pitch, dz, 17, 17, 1,
                           CylinderFilter::Orientation::Z, 0);

  Position r(-0.13232, 0.456197783422936, -32.);
  Direction u(0.43232,  -0.255, 0.8649129537704936);
  r = r - 10 * u;
  Tracker trkr(r, u, true);

  // auto index_distance = cy_filter.get_indices_tracklength(trkr, 100.);
  // std::cout << "length of the vector = " << index_distance.size() << std::endl;
  // for (auto& p : index_distance) {
  //   std::cout << "index: " << p.index[0] << ", " << p.index[0] << ", "
  //             << p.index[2] << "\tdistance = " << p.distance << std::endl;
  // }
  std::array<int, 3> on = {0, 3, 2};

  double cross_distance = 0.;
  const double sine_pol_sqr_inv = 1. / (u.x()*u.x() + u.y()*u.y());// 1./0.25192558239999996;
  auto dis_key = cy_filter.distance_to_next_index(r, u, 1./u.x(), 1./u.y(), 1./u.z(), sine_pol_sqr_inv, on, 8, 8, 0, cross_distance);
  std::cout << cross_distance << " << cross-distance" << std::endl;


  std::cout << "\n\nWe are exiting the test function. \n\n" << std::endl;
}

#endif