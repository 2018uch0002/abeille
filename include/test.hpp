#ifndef TEST_H
#define TEST_H

#include <tallies/cylinder_filter.hpp>

void test() {
  std::cout << "We are inside of the test function. \n\n" << std::endl;
  Position origin(-10.0784, -10.0784, -100.);
  const double radius = 0.4750;
  const double pitch = 1.2598;
  const double dz = 20.;

  CylinderFilter cy_filter(origin, radius, pitch, pitch, dz, 17, 17, 10,
                           CylinderFilter::Orientation::Z, 0);

  Position r(-11.4358925, 0., -154.055);
  Direction u(1E-15, 1E-15, 1.);
  Tracker trkr(r, u, true);

  auto index_distance = cy_filter.get_indices_tracklength(trkr, 100.);
  std::cout << "length of the vector = " << index_distance.size() << std::endl;
  for (auto& p : index_distance) {
    std::cout << "index: " << p.index[0] << ", " << p.index[0] << ", "
              << p.index[2] << "\tdistance = " << p.distance << std::endl;
  }
  std::array<int, 3> on = {0, 3, 2};

  double cross_distance = 10;
  cy_filter.distance_to_next_index(r, u, 1./u.x(), 1./u.y(), 1./u.z(), 2., on, 1, 2, 3, cross_distance);

  std::cout << "\n\nWe are exiting the test function. \n\n" << std::endl;
}

#endif