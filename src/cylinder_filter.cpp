#include <tallies/cylinder_filter.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

#include <cmath>
#include <sstream>
#include <vector>

CylinderFilter::CylinderFilter(Position origin, double radius, double dx,
                               double dy, double dz, std::size_t nx,
                               std::size_t ny, std::size_t nz, Orientation z_,
                               std::size_t id)
    : PositionFilter(id),
      origin_(origin),
      r_low_(),
      r_high_(),
      Nx_(nx),
      Ny_(ny),
      Nz_(nz),
      Real_nx_(nx),
      Real_ny_(ny),
      Real_nz_(nz),
      length_axis_(z_),
      radius_(radius),
      pitch_x_(dx),
      pitch_y_(dy),
      dz_(dz),
      inv_radius_(1. / radius_),
      inv_pitch_x_(),
      inv_pitch_y_(),
      inv_dz_(),
      x_index(),
      y_index(),
      z_index() {
  // Map the parameters according to the orientation
  // since this class can do the caluclation assumind z-axis as axial
  // therefore, all the parameter will be mapped to z-axis,
  // however, while return out from this class, parameters
  // will be mapped back to its original orientation
  origin_ = map_coordinate(origin);
  if (length_axis_ == Orientation::X) {
    pitch_x_ = dz;
    dz_ = dx;
    Nx_ = nz;
    Nz_ = nx;
  } else if (length_axis_ == Orientation::Y) {
    pitch_y_ = dz;
    dz_ = dy;
    Nz_ = ny;
    Ny_ = nz;
  }

  // Check for the valid radius
  if (radius_ <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " was provided with a negative or zero radius.";
    fatal_error(mssg.str());
  }

  // Make sure pitches are all >= 0.
  if (dx <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with a dx <= 0.";
    fatal_error(mssg.str());
  }

  if (dy <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with a dy <= 0.";
    fatal_error(mssg.str());
  }

  if (dz <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with a dz <= 0.";
    fatal_error(mssg.str());
  }

  // Make sure shapes are all > 0
  if (nx == 0 && length_axis_ == Orientation::X) {
    infinite_length_ = true;
  } else if (nx == 0) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with nx = 0.";
    fatal_error(mssg.str());
  }

  if (ny == 0 && length_axis_ == Orientation::Y) {
    infinite_length_ = true;
  } else if (ny == 0) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with ny = 0.";
    fatal_error(mssg.str());
  }

  if (nz == 0 && length_axis_ == Orientation::Z) {
    infinite_length_ = true;
  } else if (nz == 0) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with nz = 0.";
    fatal_error(mssg.str());
  }

  // Make sure the pitch_x and pitch_y are >= diameter
  if (pitch_x_ < 2. * radius_ || pitch_y_ < 2. * radius_) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " has dimensions which are shorter than the diameter.";
    fatal_error(mssg.str());
  }

  // Calculate inverse values
  inv_pitch_x_ = 1.0 / pitch_x_;
  inv_pitch_y_ = 1.0 / pitch_y_;
  inv_dz_ = 1.0 / dz_;

  // calculate low point for the purpose of getting the indices
  // r_low_ is mapped, means z-cooridnate will always be axial direction
  double low_x = origin_.x() - pitch_x_ * 0.5;
  double low_y = origin_.y() - pitch_y_ * 0.5;
  r_low_ = Position(low_x, low_y, origin_.z());
  r_high_ = Position(low_x + pitch_x_ * static_cast<double>(Nx_),
                     low_y + pitch_y_ * static_cast<double>(Ny_),
                     origin_.z() + dz_ * static_cast<double>(Nz_));

  // Assign the x_index, y_index, and z_index
  // these location will be based on the real_nx, real_ny, and real_nz
  x_index = 0;
  y_index = 1;
  z_index = 2;
  if (Real_nx_ == 1) {
    x_index = 0;
    y_index--;
    z_index--;
  }

  if (Real_ny_ == 1) {
    y_index = 0;
    z_index--;
  }

  if (Real_nz_ == 1) {
    z_index = 0;
  }
}

StaticVector3 CylinderFilter::get_indices(const Tracker& tktr) const {
  const Position r = map_coordinate(tktr.r());

  const int nx =
      static_cast<int>(std::floor((r.x() - r_low_.x()) * inv_pitch_x_));
  const int ny =
      static_cast<int>(std::floor((r.y() - r_low_.y()) * inv_pitch_y_));
  const int nz =
      infinite_length_
          ? 0
          : static_cast<int>(std::floor((r.z() - r_low_.z()) * inv_dz_));

  double new_origin_x = origin_.x() + pitch_x_ * static_cast<double>(nx);
  double new_origin_y = origin_.y() + pitch_y_ * static_cast<double>(ny);

  StaticVector3 indices;
  // check if nx, ny, and nz are positive
  // and less the number of bins in that direction
  if (nx >= 0 && nx < static_cast<int>(Nx_) && ny >= 0 &&
      ny < static_cast<int>(Ny_) &&
      ((nz >= 0 && nz < static_cast<int>(Nz_)) || infinite_length_)) {
    // check if the position is inside the circular radius or not
    if (std::sqrt((new_origin_x - r.x()) * (new_origin_x - r.x()) +
                  (new_origin_y - r.y()) * (new_origin_y - r.y())) <=
        (radius_ + 1E-15)) {
      indices.push_back(static_cast<std::size_t>(nx));
      indices.push_back(static_cast<std::size_t>(ny));
      indices.push_back(static_cast<std::size_t>(nz));
      map_indexes(indices);
      return reduce_dimension(indices[0], indices[1], indices[2]);
    }
  }
  return indices;
}

StaticVector3 CylinderFilter::get_position_index(const Position& r) const {
  const Position maped_r = map_coordinate(r);

  const int nx =
      static_cast<int>(std::floor((maped_r.x() - r_low_.x()) * inv_pitch_x_));
  const int ny =
      static_cast<int>(std::floor((maped_r.y() - r_low_.y()) * inv_pitch_y_));
  const int nz =
      infinite_length_
          ? 0
          : static_cast<int>(std::floor((maped_r.z() - r_low_.z()) * inv_dz_));

  double new_origin_x = origin_.x() + pitch_x_ * static_cast<double>(nx);
  double new_origin_y = origin_.y() + pitch_y_ * static_cast<double>(ny);

  StaticVector3 indices;
  // check if nx, ny, and nz are positive
  // and less the number of bins in that direction
  if (nx >= 0 && nx < static_cast<int>(Nx_) && ny >= 0 &&
      ny < static_cast<int>(Ny_) &&
      ((nz >= 0 && nz < static_cast<int>(Nz_)) || infinite_length_)) {
    // check if the position is inside the circular radius or not
    if (std::sqrt((new_origin_x - maped_r.x()) * (new_origin_x - maped_r.x()) +
                  (new_origin_y - maped_r.y()) *
                      (new_origin_y - maped_r.y())) <= (radius_ + 1E-15)) {
      indices.push_back(static_cast<std::size_t>(nx));
      indices.push_back(static_cast<std::size_t>(ny));
      indices.push_back(static_cast<std::size_t>(nz));
      map_indexes(indices);
      return reduce_dimension(indices[0], indices[1], indices[2]);
    }
  }
  return indices;
}

StaticVector3 CylinderFilter::get_shape() const {
  // shape in case of one cylinder with finite length
  if (Real_nx_ == 1 && Real_ny_ == 1 && Real_nz_ == 1) {
    return {1};
  }

  // shape in case of one infinite-cylinder
  if (infinite_length_ == true)
    if (Nx_ == 1 && Ny_ == 1) {
      return {1};
    }

  return reduce_dimension(Real_nx_, Real_ny_, Real_nz_);
}

std::vector<TracklengthDistance> CylinderFilter::get_indices_tracklength(
    const Tracker& trkr, double d_flight) const {
  if (infinite_length_) {
    fatal_error(
        "the conditions for infinite cylinder filter is not yet implemented "
        "for the track-length.");
  }
  std::vector<TracklengthDistance> indices_tracklength;
  TracklengthDistance trlen_d;

  Position r = map_coordinate(trkr.r());
  const Direction u = map_direction(trkr.u());
  const double ux_inv = 1. / u.x();
  const double uy_inv = 1. / u.y();
  const double uz_inv = 1. / u.z();

  bool inside_bin = false;

  int i = 0, j = 0, k = 0;
  std::array<int, 3> on;
  on.fill(0);  // to know we are on which tile
  initialize_indices(r, u, i, j, k, on);

  // check if particle is inside any bin.
  if ((i >= 0 && i < static_cast<int>(Nx_)) &&
      (j >= 0 && j < static_cast<int>(Ny_)) &&
      ((k >= 0 && k < static_cast<int>(Nz_)) || infinite_length_)) {
    inside_bin = true;
  } else {
    // if particle is not inside, then it can pass through the tally-region.
    if (find_entry_point(r, u, ux_inv, uy_inv, uz_inv, d_flight) == false) {
      return indices_tracklength;
    }

    initialize_indices(r, u, i, j, k, on);
    if ((i >= 0 && i < static_cast<int>(Nx_)) &&
        (j >= 0 && j < static_cast<int>(Ny_)) &&
        ((k >= 0 && k < static_cast<int>(Nz_)) || infinite_length_)) {
      inside_bin = true;
    } else {
      // This is a problem, in theory, we should now be inside the tally
      // region. We will therefore spew a warning here.
      warning("Could not locate tile after fast forward to mesh entry.\n");
    }
  }

  // Generallized method to get the distance will not work,
  // if the particle is moving perpendicular to radial plane.
  if (std::abs(1. - std::abs(u.z())) < SURFACE_COINCIDENT) {
    // since we are moving perpendicular to radial plane,
    // check first weather we are inside the circle or not.
    // if not inside the cylinder's circle, then don't socre and return
    const double new_origin_x = origin_.x() + pitch_x_ * static_cast<double>(i);
    const double new_origin_y = origin_.y() + pitch_y_ * static_cast<double>(j);
    const double xp = (r.x() - new_origin_x);
    const double yp = (r.y() - new_origin_y);

    if (xp * xp + yp * yp >= radius_ * radius_ + 1E-15) {
      std::cout << origin_.y() << "\t" << pitch_y_ << std::endl;
      std::cout << "WE COMING HERE  " << xp << "\t" << yp << "\t"
                << xp * xp + yp * yp << "\tj= " << j << "\t" << new_origin_y
                << std::endl;
      return indices_tracklength;
    }

    std::cout << "ARE WE COMING HERE \n" << std::endl;

    const int k_increment = static_cast<int>(std::copysign(1., u.z()));

    // now store the distance travelled in the first bin.
    std::size_t ui = static_cast<std::size_t>(i);
    std::size_t uj = static_cast<std::size_t>(j);
    std::size_t uk = static_cast<std::size_t>(k);
    trlen_d.index = reduce_dimension(ui, uj, uk);

    double zmin = r_low_.z() + static_cast<double>(k) * dz_;
    double cross_dist = std::abs(r.z() - zmin - (k_increment == -1 ? 0. : dz_));
    trlen_d.distance = std::min(d_flight, cross_dist);
    indices_tracklength.push_back(trlen_d);

    d_flight -= cross_dist;  // reduce flight distance by the cross distance.
    r = r + cross_dist * u;  // move the position to new position.

    if (d_flight < 0) return indices_tracklength;

    // store the distance travelled in the last bin, if it is inside bin.
    const double rz_last = r.z() + d_flight * u.z();
    const int nz_final =
        static_cast<int>(std::floor((rz_last - r_low_.z()) * inv_dz_));
    if (nz_final >= 0 && nz_final < static_cast<int>(Nz_)) {
      zmin = r_low_.z() + static_cast<double>(nz_final) * dz_;
      cross_dist = std::abs(rz_last - zmin - (k_increment == 1 ? 0 : dz_));
      trlen_d.distance = std::min(d_flight, cross_dist);

      uk = static_cast<std::size_t>(nz_final);
      trlen_d.index = reduce_dimension(ui, uj, uk);
      indices_tracklength.push_back(trlen_d);
      d_flight -= cross_dist;
    }

    // add index and distance of remaining of the scoring bins
    // so, start from the index of scoring bin after first index
    while (d_flight > 0.) {
      k += k_increment;

      if (0 <= k && k < static_cast<int>(Nz_)) {
        trlen_d.distance = std::min(d_flight, dz_);

        uk = static_cast<std::size_t>(k);
        trlen_d.index = reduce_dimension(ui, uj, uk);

        indices_tracklength.push_back(trlen_d);

      } else {
        // If we arrive here, it means that we have left the tally region.
        return indices_tracklength;
      }
      // subtract the travelled distance
      d_flight -= dz_;
      r = r + dz_ * u;
    }

    std::cout << "ARE WE COMING HERE " << indices_tracklength.size() << " \n"
              << std::endl;

    return indices_tracklength;
  }

  return indices_tracklength;
}

void CylinderFilter::initialize_indices(const Position& r, const Direction& u,
                                        int& i, int& j, int& k,
                                        std::array<int, 3>& on) const {
  on.fill(0);

  // get the index based on the position
  i = static_cast<int>(std::floor((r.x() - r_low_.x()) * inv_pitch_x_));
  j = static_cast<int>(std::floor((r.y() - r_low_.y()) * inv_pitch_y_));
  k = static_cast<int>(std::floor((r.z() - r_low_.z()) * inv_dz_));

  // Get tile boundaries
  const double xl = r_low_.x() + static_cast<double>(i) * pitch_x_;
  const double xh = xl + pitch_x_;
  const double yl = r_low_.y() + static_cast<double>(j) * pitch_y_;
  const double yh = yl + pitch_y_;
  const double zl = r_low_.z() + static_cast<double>(k) * dz_;
  const double zh = zl + dz_;

  // It is necessary to handle case of being on a tile boundary.
  if (std::abs(xl - r.x()) < SURFACE_COINCIDENT) {
    if (u.x() < 0.) {
      i--;
      on[0] = 1;
    } else {
      on[0] = -1;
    }
  } else if (std::abs(xh - r.x()) < SURFACE_COINCIDENT) {
    if (u.x() < 0.) {
      on[0] = 1;
    } else {
      i++;
      on[0] = -1;
    }
  }

  if (std::abs(yl - r.y()) < SURFACE_COINCIDENT) {
    if (u.y() < 0.) {
      j--;
      on[1] = 1;
    } else {
      on[1] = -1;
    }
  } else if (std::abs(yh - r.y()) < SURFACE_COINCIDENT) {
    if (u.y() < 0.) {
      on[1] = 1;
    } else {
      j++;
      on[1] = -1;
    }
  }

  if (std::abs(zl - r.z()) < SURFACE_COINCIDENT) {
    if (u.z() < 0.) {
      k--;
      on[2] = 1;
    } else {
      on[2] = -1;
    }
  } else if (std::abs(zh - r.z()) < SURFACE_COINCIDENT) {
    if (u.z() < 0.) {
      on[2] = 1;
    } else {
      k++;
      on[2] = -1;
    }
  }
}

bool CylinderFilter::find_entry_point(Position& r, const Direction& u,
                                      const double& ux_inv,
                                      const double& uy_inv,
                                      const double& uz_inv,
                                      double& d_flight) const {
  double d_min = (r_low_.x() - r.x()) * ux_inv;
  double d_max = (r_high_.x() - r.x()) * ux_inv;

  if (d_min > d_max) {
    std::swap(d_min, d_max);
  }

  double d_y_min = (r_low_.y() - r.y()) * uy_inv;
  double d_y_max = (r_high_.y() - r.y()) * uy_inv;

  if (d_y_min > d_y_max) {
    std::swap(d_y_min, d_y_max);
  }

  if ((d_min > d_y_max) || (d_y_min > d_max)) {
    return false;
  }

  if (d_y_min > d_min) {
    d_min = d_y_min;
  }

  if (d_y_max < d_max) {
    d_max = d_y_max;
  }

  double d_z_min = (r_low_.z() - r.z()) * uz_inv;
  double d_z_max = (r_high_.z() - r.z()) * uz_inv;

  if (d_z_min > d_z_max) {
    std::swap(d_z_min, d_z_max);
  }

  if ((d_min > d_z_max) || (d_z_min > d_max)) {
    return false;
  }

  if (d_z_min > d_min) {
    d_min = d_z_min;
  }

  if (d_z_max < d_max) {
    d_max = d_z_max;
  }

  if (d_max < d_min) {
    std::swap(d_max, d_min);
  }

  if ((d_max < 0.) && (d_min < 0.)) {
    return false;
  }

  if (d_min < 0.) {
    // If we are here, this means that r is actually inside the mesh, but is
    // really close to the edge, and we have a direction taking us out.
    // We should return false here, so that we don't score anything for this
    // particle track.
    return false;
  }

  // If we get here, we intersect the box. Let's update the position and the
  // flight distance.
  r = r + d_min * u;
  d_flight -= d_min;
  return true;
}

std::pair<double, int> CylinderFilter::distance_to_next_index(
    const Position& r, const Direction& u, const double& ux_inv,
    const double& uy_inv, const double& uz_inv, const double& sine_polar_angle,
    const std::array<int, 3>& on, int i, int j, int k, double& cross_distance) const {

    // Set our initial value for the distance and the index change
    double box_dist = INF;
    cross_distance = 0.;
    int key = 0;

    // Check all six sides
    const double diff_xl = r_low_.x() + static_cast<double>(i) * pitch_x_ - r.x();
    const double diff_xh = diff_xl + pitch_x_;
    const double diff_yl = r_low_.y() + static_cast<double>(j) * pitch_y_ - r.y();
    const double diff_yh = diff_yl + pitch_y_;
    const double diff_zl = r_low_.z() + static_cast<double>(k) * dz_ - r.z();
    const double diff_zh = diff_zl + dz_;

    const double d_xl = diff_xl * ux_inv;
    const double d_xh = diff_xh * ux_inv;
    const double d_yl = diff_yl * uy_inv;
    const double d_yh = diff_yh * uy_inv;
    const double d_zl = diff_zl * uz_inv;
    const double d_zh = diff_zh * uz_inv;

    if (d_xl > 0. && d_xl < box_dist && on[0] != -1) {
      box_dist = d_xl;
      key = -1;
    }

    if (d_xh > 0. && d_xh < box_dist && on[0] != 1) {
      box_dist = d_xh;
      key = 1;
    }

    if (d_yl > 0. && d_yl < box_dist && on[1] != -1) {
      box_dist = d_yl;
      key = -2;
    }

    if (d_yh > 0. && d_yh < box_dist && on[1] != 1) {
      box_dist = d_yh;
      key = 2;
    }

    if (d_zl > 0. && d_zl < box_dist && on[2] != -1) {
      box_dist = d_zl;
      key = -3;
    }

    if (d_zh > 0. && d_zh < box_dist && on[2] != 1) {
      box_dist = d_zh;
      key = 3;
    }

    return {box_dist, key};
  }

double CylinderFilter::z_min(const StaticVector3& index) const {
  // note that "index" is not orientated according to class in general
  if (Real_nz_ == 1) {
    return r_low_.z();
  }
  return (r_low_.z() + static_cast<double>(index[z_index]) * dz_);
}

double CylinderFilter::z_max(const StaticVector3& index) const {
  // note that "index" is not orientated according to class in general
  if (Real_nz_ == 1) {
    return r_low_.z() + dz_;
  }
  return (r_low_.z() + static_cast<double>(index[z_index]) * dz_ + dz_);
}

Position CylinderFilter::get_center(const StaticVector3& index,
                                    const bool is_map = true) const {
  // note that "index" is not orientated according to class in general
  double new_origin_x = origin_.x();
  double new_origin_y = origin_.y();
  double new_origin_z = origin_.z();
  if (Real_nx_ > 1) {
    new_origin_x += pitch_x_ * static_cast<double>(index[x_index]);
  }

  if (Real_ny_ > 1) {
    new_origin_y += pitch_y_ * static_cast<double>(index[y_index]);
  }

  if (Real_nz_ > 1) {
    new_origin_z += dz_ * static_cast<double>(index[z_index]);
  }

  // if is_map is false, the send the Position based on the class-orientation
  if (is_map == false) {
    return Position(new_origin_x, new_origin_y, new_origin_z);
  }
  return map_coordinate(Position(new_origin_x, new_origin_y, new_origin_z));
}

// first will be the scaled radius and second will be the angle
std::pair<double, double> CylinderFilter::get_scaled_radius_and_angle(
    const StaticVector3& indices, const Position& r) const {
  // the Position is not according to the class, so map it.
  Position mapped_r = map_coordinate(r);
  // get the new-origin cooredinates based on the indices
  // the new-origin should be according to the class orientation
  Position new_origin = get_center(indices, false);

  const double base = mapped_r.x() - new_origin.x();
  const double height = mapped_r.y() - new_origin.y();
  // scaled-radius
  const double scaled_r =
      (std::sqrt(base * base + height * height)) * inv_radius_;
  // theta
  double theta = std::atan2(height, base);
  // atan2 provides the anlges between [-PI, PI], instead of [0, 2*PI]
  if (theta < 0.) {
    return {scaled_r, theta + 2 * PI};
  }
  return {scaled_r, theta};
}

void CylinderFilter::write_to_hdf5(H5::Group& grp) const {
  // Save id in attributes
  if (grp.hasAttribute("id")) {
    grp.deleteAttribute("id");
  }
  grp.createAttribute("id", this->id());

  // Save type in attributes
  if (grp.hasAttribute("type")) {
    grp.deleteAttribute("type");
  }
  grp.createAttribute("type", "cylinder-filter");

  // Save origin position
  std::array<double, 3> origin{origin_.x(), origin_.y(), origin_.z()};
  if (grp.hasAttribute("origin")) {
    grp.deleteAttribute("origin");
  }
  grp.createAttribute("origin", origin);

  // Save radius
  if (grp.hasAttribute("radius")) {
    grp.deleteAttribute("radius");
  }
  grp.createAttribute("radius", this->radius());

  // Save axis
  if (grp.hasAttribute("axis")) {
    grp.deleteAttribute("axis");
  }
  switch (length_axis_) {
    case Orientation::X:
      grp.createAttribute("axis", "x");
      break;

    case Orientation::Y:
      grp.createAttribute("axis", "y");
      break;

    case Orientation::Z:
      grp.createAttribute("axis", "z");
      break;
  }

  // Save shape
  std::array<std::size_t, 3> shape{Real_nx_, Real_ny_, Real_nz_};
  if (grp.hasAttribute("shape")) {
    grp.deleteAttribute("shape");
  }
  grp.createAttribute("shape", shape);

  // Save pitch
  std::array<double, 3> pitch{pitch_x_, pitch_y_, dz_};
  if (length_axis_ == Orientation::X) {
    pitch[0] = dz_;
    pitch[2] = pitch_x_;
  } else if (length_axis_ == Orientation::Y) {
    pitch[1] = dz_;
    pitch[2] = pitch_y_;
  }
  if (grp.hasAttribute("pitch")) {
    grp.deleteAttribute("pitch");
  }
  grp.createAttribute("pitch", pitch);
}

std::shared_ptr<CylinderFilter> make_cylinder_filter(const YAML::Node& node) {
  // Get the id
  if (!node["id"] || node["id"].IsScalar() == false) {
    fatal_error("Invalid id is given for the position-filter.");
  }
  std::size_t id = node["id"].as<std::size_t>();

  // check and get the origin
  if (!node["origin"] || node["origin"].IsSequence() == false ||
      node["origin"].size() != 3) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " is missing a valid origin entry.";
    fatal_error(mssg.str());
  }
  std::vector<double> origin_point = node["origin"].as<std::vector<double>>();
  Position origin(origin_point[0], origin_point[1], origin_point[2]);

  // Get the radius
  if (!node["radius"] || node["radius"].IsScalar() == false) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " is missing a valid radius entry.";
    fatal_error(mssg.str());
  }
  const double radius = node["radius"].as<double>();
  if (radius <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " has a radius which is <= 0.";
    fatal_error(mssg.str());
  }

  // Get the axial-direction axis name
  if (!node["axis"] || node["axis"].IsScalar() == false) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " is missing a valid axis entry.";
    fatal_error(mssg.str());
  }
  const std::string axial_axis = node["axis"].as<std::string>();
  if (axial_axis != "x" && axial_axis != "y" && axial_axis != "z") {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " was provided with invalid axis entry \"" << axial_axis << "\".";
    fatal_error(mssg.str());
  }
  CylinderFilter::Orientation orientation;
  if (axial_axis == "x")
    orientation = CylinderFilter::Orientation::X;
  else if (axial_axis == "y")
    orientation = CylinderFilter::Orientation::Y;
  else
    orientation = CylinderFilter::Orientation::Z;

  // The default shape will nx = 1, ny = 1, and nz = 1;
  std::size_t nx = 1;
  std::size_t ny = 1;
  std::size_t nz = 1;
  // Get the shape
  if (node["shape"] && node["shape"].IsSequence() &&
      node["shape"].size() == 3) {
    std::vector<std::size_t> shape =
        node["shape"].as<std::vector<std::size_t>>();
    nx = shape[0];
    ny = shape[1];
    nz = shape[2];
  } else if (node["shape"]) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " has an invalid shape entry.";
    fatal_error(mssg.str());
  }

  // Get the pitch of the lattice/box.
  if (!node["pitch"] || node["pitch"].IsSequence() == false ||
      node["pitch"].size() != 3) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " has invalid pitch entry.";
    fatal_error(mssg.str());
  }
  std::vector<double> pitches = node["pitch"].as<std::vector<double>>();

  double pitch_x, pitch_y, length;
  if (axial_axis == "z") {
    pitch_x = pitches[0];
    pitch_y = pitches[1];
    length = pitches[2];
  } else if (axial_axis == "x") {
    pitch_x = pitches[2];
    pitch_y = pitches[1];
    length = pitches[0];
  } else {
    pitch_x = pitches[0];
    pitch_y = pitches[2];
    length = pitches[1];
  }

  if (pitch_x <= 0. || pitch_y <= 0. || length <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " has pitches which are <= 0.";
    fatal_error(mssg.str());
  }

  if (pitch_x < 2. * radius || pitch_y < 2. * radius) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " has pitches which are too small for the given radius.";
    fatal_error(mssg.str());
  }

  // Make the filter
  return std::make_shared<CylinderFilter>(origin, radius, pitches[0],
                                          pitches[1], pitches[2], nx, ny, nz,
                                          orientation, id);
}