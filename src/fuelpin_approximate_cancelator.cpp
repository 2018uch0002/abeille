#include <cancelator/fuelpin_approximate_cancelator.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>

#include <xtensor/xarray.hpp>

#include <algorithm>
#include <cmath>

// If the number of rings are given.
FuelPinApproxCancelator::FuelPinApproxCancelator(
    Position low, Position high, std::vector<double> assembly_pitch,
    std::vector<std::uint32_t> assembly_shape,
    std::vector<double> fuelpin_pitch,
    std::vector<std::uint32_t> fuelpins_per_assembly,
    std::vector<double> inter_assembly_gap, double radius, std::uint32_t Nr,
    Orientation axial_orientation, bool loop)
    : bins_(),
      energy_edges_(),
      r_low_(low),
      r_high_(high),
      assembly_pitch_x_(),
      assembly_pitch_y_(),
      assembly_pitch_z_(),
      assembly_pitch_x_inv_(),
      assembly_pitch_y_inv_(),
      assembly_pitch_z_inv_(),
      pitch_x_(),
      pitch_y_(),
      pitch_z_(),
      pitch_x_inv_(),
      pitch_y_inv_(),
      pitch_z_inv_(),
      inter_asmbly_gap_x_(),
      inter_asmbly_gap_y_(),
      inter_asmbly_gap_z_(),
      radius_ring_square_(),
      Nx_asmbly_(),
      Ny_asmbly_(),
      Nz_asmbly_(),
      Nx_per_asmbly_(),
      Ny_per_asmbly_(),
      Nz_per_asmbly_(),
      Nr_(Nr),
      Ns_(4),
      Ne_(1),
      length_axis_(axial_orientation),
      loop(loop) {
  // check axial orientation
  if (length_axis_ != Orientation::Z) {
    fatal_error(
        "Currently only fuelpin with axial direction along the z-direction is "
        "supported.");
  }
  // check low-point is lower than the high-point
  if (r_low_.x() >= r_high_.x() || r_low_.y() >= r_high_.y() ||
      r_low_.z() >= r_high_.z()) {
    fatal_error(
        "Low position is not lower than high position in "
        "FuelPinApproxCancelator.");
  }

  // check the assembly-pitch
  if (assembly_pitch.size() != 3) {
    fatal_error("Size of the assembly-pitch must be 3.");
  }
  assembly_pitch_x_ = assembly_pitch[0];
  assembly_pitch_y_ = assembly_pitch[1];
  assembly_pitch_z_ = assembly_pitch[2];
  if (assembly_pitch_x_ <= 0. || assembly_pitch_y_ <= 0. ||
      assembly_pitch_z_ <= 0.) {
    fatal_error("Assembly-pitch must be positive and > 0 in each direction.");
  }

  assembly_pitch_x_inv_ = 1. / assembly_pitch_x_;
  assembly_pitch_y_inv_ = 1. / assembly_pitch_y_;
  assembly_pitch_z_inv_ = 1. / assembly_pitch_z_;

  // check the assembly shape
  if (assembly_shape.size() != 3) {
    fatal_error("Size of the assembly-shape must be 3.");
  }
  Nx_asmbly_ = assembly_shape[0];
  Ny_asmbly_ = assembly_shape[1];
  Nz_asmbly_ = assembly_shape[2];

  // check the fuel-pin-pitch
  if (fuelpin_pitch.size() != 3) {
    fatal_error("Size of the fuel_pin_pitch must be 3.");
  }

  pitch_x_ = fuelpin_pitch[0];
  pitch_y_ = fuelpin_pitch[1];
  pitch_z_ = fuelpin_pitch[2];
  if (pitch_x_ <= 0. || pitch_y_ <= 0. || pitch_z_ <= 0.) {
    fatal_error("Fuelpin-pitch must be positive and > 0 in each direction.");
  }

  pitch_x_inv_ = 1. / pitch_x_;
  pitch_y_inv_ = 1. / pitch_y_;
  pitch_z_inv_ = 1. / pitch_z_;

  // check the fuelpin shape within a assembly
  if (fuelpins_per_assembly.size() != 3) {
    fatal_error("Size of the shape of fuelpins per assembly must be 3.");
  }
  Nx_per_asmbly_ = fuelpins_per_assembly[0];
  Ny_per_asmbly_ = fuelpins_per_assembly[1];
  Nz_per_asmbly_ = fuelpins_per_assembly[2];

  // get the inter-assembly-gaps
  if (inter_assembly_gap.size() != 3) {
    fatal_error("Size of the inter-assebly-gap must be 3.");
  }
  inter_asmbly_gap_x_ = inter_assembly_gap[0];
  inter_asmbly_gap_y_ = inter_assembly_gap[1];
  inter_asmbly_gap_z_ = inter_assembly_gap[2];
  if (inter_asmbly_gap_x_ < 0. || inter_asmbly_gap_y_ < 0. ||
      inter_asmbly_gap_z_ < 0.) {
    fatal_error("Inter-assebly-gap must be >= 0 in each direction.");
  }

  // get the equi-arial (or volume) rings
  radius_ring_square_.reserve(Nr_);
  double r1 = radius / std::sqrt(static_cast<double>(Nr_));
  radius_ring_square_.push_back(r1 * r1);

  for (std::size_t ir = 1; ir < Nr_; ir++) {
    const double r_ir = r1 * std::sqrt(static_cast<double>(ir + 1));
    radius_ring_square_.push_back(r_ir * r_ir);
  }

  // set the strides
  const std::uint32_t Ny = Ny_asmbly_ * Ny_per_asmbly_;
  const std::uint32_t Nz = Nz_asmbly_ * Nz_per_asmbly_;
  Si_ = Ne_ * Ns_ * Nr_ * Nz * Ny;
  Sj_ = Ne_ * Ns_ * Nr_ * Nz;
  Sk_ = Ne_ * Ns_ * Nr_;
  Sr_ = Ne_ * Ns_;
  Ss_ = Ne_;
  Se_ = 1;
}

// If the number of rings and the energy-groups are given.
FuelPinApproxCancelator::FuelPinApproxCancelator(
    Position low, Position high, std::vector<double> assembly_pitch,
    std::vector<std::uint32_t> assembly_shape,
    std::vector<double> fuelpin_pitch,
    std::vector<std::uint32_t> fuelpins_per_assembly,
    std::vector<double> inter_assembly_gap, double radius, std::uint32_t Nr,
    Orientation axial_orientation, std::vector<double> energy_bounds, bool loop)
    : bins_(),
      energy_edges_(energy_bounds),
      r_low_(low),
      r_high_(high),
      assembly_pitch_x_(),
      assembly_pitch_y_(),
      assembly_pitch_z_(),
      assembly_pitch_x_inv_(),
      assembly_pitch_y_inv_(),
      assembly_pitch_z_inv_(),
      pitch_x_(),
      pitch_y_(),
      pitch_z_(),
      pitch_x_inv_(),
      pitch_y_inv_(),
      pitch_z_inv_(),
      inter_asmbly_gap_x_(),
      inter_asmbly_gap_y_(),
      inter_asmbly_gap_z_(),
      radius_ring_square_(),
      Nx_asmbly_(),
      Ny_asmbly_(),
      Nz_asmbly_(),
      Nx_per_asmbly_(),
      Ny_per_asmbly_(),
      Nz_per_asmbly_(),
      Nr_(Nr),
      Ns_(4),
      Ne_(1),
      length_axis_(axial_orientation),
      loop(loop) {
  // check axial orientation
  if (length_axis_ != Orientation::Z) {
    fatal_error(
        "Currently only fuelpin with axial direction along the z-direction is "
        "supported.");
  }
  // check low-point is lower than the high-point
  if (r_low_.x() >= r_high_.x() || r_low_.y() >= r_high_.y() ||
      r_low_.z() >= r_high_.z()) {
    fatal_error(
        "Low position is not lower than high position in "
        "FuelPinApproxCancelator.");
  }

  // Make sure energy bins are valid
  if (energy_edges_.size() < 2) {
    fatal_error(
        "energy-edges must have at least two entries in "
        "FuelPinApproxCancelator.");
  }

  if (!std::is_sorted(energy_edges_.begin(), energy_edges_.end())) {
    fatal_error("energy-edges must be sorted in FuelPinApproxCancelator.");
  }

  if (energy_edges_.front() < 0.) {
    fatal_error(
        "All energy-edges must be greater than or equal to zero in "
        "FuelPinApproxCancelator.");
  }

  // check the assembly-pitch
  if (assembly_pitch.size() != 3) {
    fatal_error("Size of the assembly-pitch must be 3.");
  }
  assembly_pitch_x_ = assembly_pitch[0];
  assembly_pitch_y_ = assembly_pitch[1];
  assembly_pitch_z_ = assembly_pitch[2];
  if (assembly_pitch_x_ <= 0. || assembly_pitch_y_ <= 0. ||
      assembly_pitch_z_ <= 0.) {
    fatal_error("Assembly-pitch must be positive and > 0 in each direction.");
  }

  assembly_pitch_x_inv_ = 1. / assembly_pitch_x_;
  assembly_pitch_y_inv_ = 1. / assembly_pitch_y_;
  assembly_pitch_z_inv_ = 1. / assembly_pitch_z_;

  // check the assembly shape
  if (assembly_shape.size() != 3) {
    fatal_error("Size of the assembly-shape must be 3.");
  }
  Nx_asmbly_ = assembly_shape[0];
  Ny_asmbly_ = assembly_shape[1];
  Nz_asmbly_ = assembly_shape[2];

  // check the fuel-pin-pitch
  if (fuelpin_pitch.size() != 3) {
    fatal_error("Size of the fuel_pin_pitch must be 3.");
  }

  pitch_x_ = fuelpin_pitch[0];
  pitch_y_ = fuelpin_pitch[1];
  pitch_z_ = fuelpin_pitch[2];
  if (pitch_x_ <= 0. || pitch_y_ <= 0. || pitch_z_ <= 0.) {
    fatal_error("Fuelpin-pitch must be positive and > 0 in each direction.");
  }

  pitch_x_inv_ = 1. / pitch_x_;
  pitch_y_inv_ = 1. / pitch_y_;
  pitch_z_inv_ = 1. / pitch_z_;

  // check the fuelpin shape within a assembly
  if (fuelpins_per_assembly.size() != 3) {
    fatal_error("Size of the shape of fuelpins per assembly must be 3.");
  }
  Nx_per_asmbly_ = fuelpins_per_assembly[0];
  Ny_per_asmbly_ = fuelpins_per_assembly[1];
  Nz_per_asmbly_ = fuelpins_per_assembly[2];

  // get the inter-assembly-gaps
  if (inter_assembly_gap.size() != 3) {
    fatal_error("Size of the inter-assebly-gap must be 3.");
  }
  inter_asmbly_gap_x_ = inter_assembly_gap[0];
  inter_asmbly_gap_y_ = inter_assembly_gap[1];
  inter_asmbly_gap_z_ = inter_assembly_gap[2];
  if (inter_asmbly_gap_x_ < 0. || inter_asmbly_gap_y_ < 0. ||
      inter_asmbly_gap_z_ < 0.) {
    fatal_error("Inter-assebly-gap must be >= 0 in each direction.");
  }

  // get the equi-arial (or volume) rings
  radius_ring_square_.reserve(Nr_);
  double r1 = radius / std::sqrt(static_cast<double>(Nr_));
  radius_ring_square_.push_back(r1 * r1);

  for (std::size_t ir = 1; ir < Nr_; ir++) {
    const double r_ir = r1 * std::sqrt(static_cast<double>(ir + 1));
    radius_ring_square_.push_back(r_ir * r_ir);
  }

  // set the strides
  const std::uint32_t Ny = Ny_asmbly_ * Ny_per_asmbly_;
  const std::uint32_t Nz = Nz_asmbly_ * Nz_per_asmbly_;
  Si_ = Ne_ * Ns_ * Nr_ * Nz * Ny;
  Sj_ = Ne_ * Ns_ * Nr_ * Nz;
  Sk_ = Ne_ * Ns_ * Nr_;
  Sr_ = Ne_ * Ns_;
  Ss_ = Ne_;
  Se_ = 1;
}

// If the radius of the rings are given.
FuelPinApproxCancelator::FuelPinApproxCancelator(
    Position low, Position high, std::vector<double> assembly_pitch,
    std::vector<std::uint32_t> assembly_shape,
    std::vector<double> fuelpin_pitch,
    std::vector<std::uint32_t> fuelpins_per_assembly,
    std::vector<double> inter_assembly_gap, std::vector<double> radii,
    Orientation axial_orientation, bool loop)

    : bins_(),
      energy_edges_(),
      r_low_(low),
      r_high_(high),
      assembly_pitch_x_(),
      assembly_pitch_y_(),
      assembly_pitch_z_(),
      assembly_pitch_x_inv_(),
      assembly_pitch_y_inv_(),
      assembly_pitch_z_inv_(),
      pitch_x_(),
      pitch_y_(),
      pitch_z_(),
      pitch_x_inv_(),
      pitch_y_inv_(),
      pitch_z_inv_(),
      inter_asmbly_gap_x_(),
      inter_asmbly_gap_y_(),
      inter_asmbly_gap_z_(),
      radius_ring_square_(radii),
      Nx_asmbly_(),
      Ny_asmbly_(),
      Nz_asmbly_(),
      Nx_per_asmbly_(),
      Ny_per_asmbly_(),
      Nz_per_asmbly_(),
      Nr_(),
      Ns_(4),
      Ne_(1),
      length_axis_(axial_orientation),
      loop(loop) {
  // check axial orientation
  if (length_axis_ != Orientation::Z) {
    fatal_error(
        "Currently only fuelpin with axial direction along the z-direction is "
        "supported.");
  }
  // check low-point is lower than the high-point
  if (r_low_.x() >= r_high_.x() || r_low_.y() >= r_high_.y() ||
      r_low_.z() >= r_high_.z()) {
    fatal_error(
        "Low position is not lower than high position in "
        "FuelPinApproxCancelator.");
  }

  // check the assembly-pitch
  if (assembly_pitch.size() != 3) {
    fatal_error("Size of the assembly-pitch must be 3.");
  }
  assembly_pitch_x_ = assembly_pitch[0];
  assembly_pitch_y_ = assembly_pitch[1];
  assembly_pitch_z_ = assembly_pitch[2];
  if (assembly_pitch_x_ <= 0. || assembly_pitch_y_ <= 0. ||
      assembly_pitch_z_ <= 0.) {
    fatal_error("Assembly-pitch must be positive and > 0 in each direction.");
  }

  assembly_pitch_x_inv_ = 1. / assembly_pitch_x_;
  assembly_pitch_y_inv_ = 1. / assembly_pitch_y_;
  assembly_pitch_z_inv_ = 1. / assembly_pitch_z_;

  // check the assembly shape
  if (assembly_shape.size() != 3) {
    fatal_error("Size of the assembly-shape must be 3.");
  }
  Nx_asmbly_ = assembly_shape[0];
  Ny_asmbly_ = assembly_shape[1];
  Nz_asmbly_ = assembly_shape[2];

  // check the fuel-pin-pitch
  if (fuelpin_pitch.size() != 3) {
    fatal_error("Size of the fuel_pin_pitch must be 3.");
  }

  pitch_x_ = fuelpin_pitch[0];
  pitch_y_ = fuelpin_pitch[1];
  pitch_z_ = fuelpin_pitch[2];
  if (pitch_x_ <= 0. || pitch_y_ <= 0. || pitch_z_ <= 0.) {
    fatal_error("Fuelpin-pitch must be positive and > 0 in each direction.");
  }

  pitch_x_inv_ = 1. / pitch_x_;
  pitch_y_inv_ = 1. / pitch_y_;
  pitch_z_inv_ = 1. / pitch_z_;

  // check the fuelpin shape within a assembly
  if (fuelpins_per_assembly.size() != 3) {
    fatal_error("Size of the shape of fuelpins per assembly must be 3.");
  }
  Nx_per_asmbly_ = fuelpins_per_assembly[0];
  Ny_per_asmbly_ = fuelpins_per_assembly[1];
  Nz_per_asmbly_ = fuelpins_per_assembly[2];

  // get the inter-assembly-gaps
  if (inter_assembly_gap.size() != 3) {
    fatal_error("Size of the inter-assebly-gap must be 3.");
  }
  inter_asmbly_gap_x_ = inter_assembly_gap[0];
  inter_asmbly_gap_y_ = inter_assembly_gap[1];
  inter_asmbly_gap_z_ = inter_assembly_gap[2];
  if (inter_asmbly_gap_x_ < 0. || inter_asmbly_gap_y_ < 0. ||
      inter_asmbly_gap_z_ < 0.) {
    fatal_error("Inter-assebly-gap must be >= 0 in each direction.");
  }

  // get the radii of the different rings
  Nr_ = static_cast<std::uint32_t>(radius_ring_square_.size());
  if (Nr_ == 0) {
    fatal_error("Size of radii must greater than 0.");
  }
  if (!std::is_sorted(radius_ring_square_.begin(), radius_ring_square_.end())) {
    fatal_error("radius_ring must be sorted in FuelPinApproxCancelator.");
  }

  if (radius_ring_square_.front() < 0.) {
    fatal_error(
        "All radius_ring_ must be greater than or equal to zero in "
        "FuelPinApproxCancelator.");
  }

  // store the square of the give radius.
  for (std::size_t ir = 0; ir < Nr_; ir++) {
    const double rad = radius_ring_square_[ir];
    radius_ring_square_[ir] = rad * rad;
  }

  // set the strides
  const std::uint32_t Ny = Ny_asmbly_ * Ny_per_asmbly_;
  const std::uint32_t Nz = Nz_asmbly_ * Nz_per_asmbly_;
  Si_ = Ne_ * Ns_ * Nr_ * Nz * Ny;
  Sj_ = Ne_ * Ns_ * Nr_ * Nz;
  Sk_ = Ne_ * Ns_ * Nr_;
  Sr_ = Ne_ * Ns_;
  Ss_ = Ne_;
  Se_ = 1;
}

// If the radius of the rings and the energy-groups are given.
FuelPinApproxCancelator::FuelPinApproxCancelator(
    Position low, Position high, std::vector<double> assembly_pitch,
    std::vector<std::uint32_t> assembly_shape,
    std::vector<double> fuelpin_pitch,
    std::vector<std::uint32_t> fuelpins_per_assembly,
    std::vector<double> inter_assembly_gap, std::vector<double> radii,
    Orientation axial_orientation, std::vector<double> energy_bounds, bool loop)

    : bins_(),
      energy_edges_(energy_bounds),
      r_low_(low),
      r_high_(high),
      assembly_pitch_x_(),
      assembly_pitch_y_(),
      assembly_pitch_z_(),
      assembly_pitch_x_inv_(),
      assembly_pitch_y_inv_(),
      assembly_pitch_z_inv_(),
      pitch_x_(),
      pitch_y_(),
      pitch_z_(),
      pitch_x_inv_(),
      pitch_y_inv_(),
      pitch_z_inv_(),
      inter_asmbly_gap_x_(),
      inter_asmbly_gap_y_(),
      inter_asmbly_gap_z_(),
      radius_ring_square_(radii),
      Nx_asmbly_(),
      Ny_asmbly_(),
      Nz_asmbly_(),
      Nx_per_asmbly_(),
      Ny_per_asmbly_(),
      Nz_per_asmbly_(),
      Nr_(),
      Ns_(4),
      Ne_(1),
      length_axis_(axial_orientation),
      loop(loop) {
  // check axial orientation
  if (length_axis_ != Orientation::Z) {
    fatal_error(
        "Currently only fuelpin with axial direction along the z-direction is "
        "supported.");
  }
  // check low-point is lower than the high-point
  if (r_low_.x() >= r_high_.x() || r_low_.y() >= r_high_.y() ||
      r_low_.z() >= r_high_.z()) {
    fatal_error(
        "Low position is not lower than high position in "
        "FuelPinApproxCancelator.");
  }

  // Make sure energy bins are valid
  if (energy_edges_.size() < 2) {
    fatal_error(
        "energy-edges must have at least two entries in "
        "FuelPinApproxCancelator.");
  }

  if (!std::is_sorted(energy_edges_.begin(), energy_edges_.end())) {
    fatal_error("energy-edges must be sorted in FuelPinApproxCancelator.");
  }

  if (energy_edges_.front() < 0.) {
    fatal_error(
        "All energy-edges must be greater than or equal to zero in "
        "FuelPinApproxCancelator.");
  }

  // check the assembly-pitch
  if (assembly_pitch.size() != 3) {
    fatal_error("Size of the assembly-pitch must be 3.");
  }
  assembly_pitch_x_ = assembly_pitch[0];
  assembly_pitch_y_ = assembly_pitch[1];
  assembly_pitch_z_ = assembly_pitch[2];
  if (assembly_pitch_x_ <= 0. || assembly_pitch_y_ <= 0. ||
      assembly_pitch_z_ <= 0.) {
    fatal_error("Assembly-pitch must be positive and > 0 in each direction.");
  }

  assembly_pitch_x_inv_ = 1. / assembly_pitch_x_;
  assembly_pitch_y_inv_ = 1. / assembly_pitch_y_;
  assembly_pitch_z_inv_ = 1. / assembly_pitch_z_;

  // check the assembly shape
  if (assembly_shape.size() != 3) {
    fatal_error("Size of the assembly-shape must be 3.");
  }
  Nx_asmbly_ = assembly_shape[0];
  Ny_asmbly_ = assembly_shape[1];
  Nz_asmbly_ = assembly_shape[2];

  // check the fuel-pin-pitch
  if (fuelpin_pitch.size() != 3) {
    fatal_error("Size of the fuel_pin_pitch must be 3.");
  }

  pitch_x_ = fuelpin_pitch[0];
  pitch_y_ = fuelpin_pitch[1];
  pitch_z_ = fuelpin_pitch[2];
  if (pitch_x_ <= 0. || pitch_y_ <= 0. || pitch_z_ <= 0.) {
    fatal_error("Fuelpin-pitch must be positive and > 0 in each direction.");
  }

  pitch_x_inv_ = 1. / pitch_x_;
  pitch_y_inv_ = 1. / pitch_y_;
  pitch_z_inv_ = 1. / pitch_z_;

  // check the fuelpin shape within a assembly
  if (fuelpins_per_assembly.size() != 3) {
    fatal_error("Size of the shape of fuelpins per assembly must be 3.");
  }
  Nx_per_asmbly_ = fuelpins_per_assembly[0];
  Ny_per_asmbly_ = fuelpins_per_assembly[1];
  Nz_per_asmbly_ = fuelpins_per_assembly[2];

  // get the inter-assembly-gaps
  if (inter_assembly_gap.size() != 3) {
    fatal_error("Size of the inter-assebly-gap must be 3.");
  }
  inter_asmbly_gap_x_ = inter_assembly_gap[0];
  inter_asmbly_gap_y_ = inter_assembly_gap[1];
  inter_asmbly_gap_z_ = inter_assembly_gap[2];
  if (inter_asmbly_gap_x_ < 0. || inter_asmbly_gap_y_ < 0. ||
      inter_asmbly_gap_z_ < 0.) {
    fatal_error("Inter-assebly-gap must be >= 0 in each direction.");
  }

  // get the radii of the different rings
  Nr_ = static_cast<std::uint32_t>(radius_ring_square_.size());
  if (Nr_ == 0) {
    fatal_error("Size of radii must greater than 0.");
  }
  if (!std::is_sorted(radius_ring_square_.begin(), radius_ring_square_.end())) {
    fatal_error("radius_ring must be sorted in FuelPinApproxCancelator.");
  }

  if (radius_ring_square_.front() < 0.) {
    fatal_error(
        "All radius_ring_ must be greater than or equal to zero in "
        "FuelPinApproxCancelator.");
  }

  // store the square of the give radius.
  for (std::size_t ir = 0; ir < Nr_; ir++) {
    const double rad = radius_ring_square_[ir];
    radius_ring_square_[ir] = rad * rad;
  }

  // set the strides
  const std::uint32_t Ny = Ny_asmbly_ * Ny_per_asmbly_;
  const std::uint32_t Nz = Nz_asmbly_ * Nz_per_asmbly_;
  Si_ = Ne_ * Ns_ * Nr_ * Nz * Ny;
  Sj_ = Ne_ * Ns_ * Nr_ * Nz;
  Sk_ = Ne_ * Ns_ * Nr_;
  Sr_ = Ne_ * Ns_;
  Ss_ = Ne_;
  Se_ = 1;
}

void FuelPinApproxCancelator::write_output_info(H5::Group& grp) const {
  grp.createAttribute("type", "fuelpin-approximate");

  const std::array<double, 3> r_low{r_low_.x(), r_low_.y(), r_low_.z()};
  grp.createAttribute("low", r_low);

  const std::array<double, 3> r_high{r_high_.x(), r_high_.y(), r_high_.z()};
  grp.createAttribute("high", r_high);

  const std::array<double, 3> assembly_pitch{
      assembly_pitch_x_, assembly_pitch_y_, assembly_pitch_z_};
  grp.createAttribute("assembly-pitch", assembly_pitch);

  const std::array<std::size_t, 3> assembly_shape{Nx_asmbly_, Ny_asmbly_,
                                                  Nz_asmbly_};
  grp.createAttribute("assembly-shape", assembly_shape);

  const std::array<double, 3> fuelpin_pitch{pitch_x_, pitch_y_, pitch_z_};
  grp.createAttribute("fuelpin-pitch", fuelpin_pitch);

  const std::array<std::size_t, 3> fuelpin_per_assembly{
      Nx_per_asmbly_, Ny_per_asmbly_, Nz_per_asmbly_};
  grp.createAttribute("fuelpin-per-assembly", fuelpin_per_assembly);

  std::vector<double> radii;
  radii.reserve(radius_ring_square_.size());
  for (std::size_t ir = 0; ir < Nr_; ir++) {
    const double rad_sqr = radius_ring_square_[ir];
    radii.push_back(std::sqrt(rad_sqr));
  }
  grp.createAttribute("radii", radii);

  grp.createAttribute("angular-segment", Ns_);

  if (length_axis_ == Orientation::X) {
    grp.createAttribute("axial-orientation", "x");
  } else if (length_axis_ == Orientation::Y) {
    grp.createAttribute("axial-orientation", "y");
  } else if (length_axis_ == Orientation::Z) {
    grp.createAttribute("axial-orientation", "z");
  }

  if (energy_edges_.empty() == false) {
    grp.createAttribute("energy-bounds", energy_edges_);
  }
}

bool FuelPinApproxCancelator::add_particle(BankedParticle& p) {
  // Get energy index with linear search
  int l = -1;
  if (energy_edges_.empty() == false) {
    for (size_t e = 0; e < energy_edges_.size() - 1; e++) {
      if (energy_edges_[e] <= p.E && p.E <= energy_edges_[e + 1]) {
        l = static_cast<int>(e);
        break;
      }
    }
  } else {
    l = 0;
  }

  // if the energy-index is out of bounds, then don't keep the particle
  if (l < 0 || static_cast<int>(Ne_) <= l) {
    return false;
  }

  // Find the indices of the assembly based on the particle location
  int i_asmbly = static_cast<int>(
      std::floor((p.r.x() - r_low_.x()) * assembly_pitch_x_inv_));
  int j_asmbly = static_cast<int>(
      std::floor((p.r.y() - r_low_.y()) * assembly_pitch_y_inv_));
  int k_asmbly = static_cast<int>(
      std::floor((p.r.z() - r_low_.z()) * assembly_pitch_z_inv_));

  // if the assembly index is out of bounds, then don't keep the particle
  if (i_asmbly < 0 || static_cast<int>(Nx_asmbly_) <= i_asmbly ||
      j_asmbly < 0 || static_cast<int>(Ny_asmbly_) <= j_asmbly ||
      k_asmbly < 0 || static_cast<int>(Nz_asmbly_) <= k_asmbly) {
    return false;
  }

  // once it is known that in which assembly particle is present, then
  // get the location of the fuel pin
  const double assembly_x0 = r_low_.x() +
                             assembly_pitch_x_ * static_cast<double>(i_asmbly) +
                             0.5 * inter_asmbly_gap_x_;
  const double assembly_y0 = r_low_.y() +
                             assembly_pitch_y_ * static_cast<double>(j_asmbly) +
                             0.5 * inter_asmbly_gap_y_;
  const double assembly_z0 = r_low_.z() +
                             assembly_pitch_z_ * static_cast<double>(k_asmbly) +
                             0.5 * inter_asmbly_gap_z_;

  // get the index of the fuel-pin within that assembly
  int i = static_cast<int>(std::floor((p.r.x() - assembly_x0) * pitch_x_inv_));
  int j = static_cast<int>(std::floor((p.r.y() - assembly_y0) * pitch_y_inv_));
  int k = static_cast<int>(std::floor((p.r.z() - assembly_z0) * pitch_z_inv_));

  // if the particle are happens to be inside the inter assembly gap, then
  // the index will be out of bounds, so don't keep the particle
  if (i < 0 || static_cast<int>(Nx_per_asmbly_) <= i || j < 0 ||
      static_cast<int>(Ny_per_asmbly_) <= j || k < 0 ||
      static_cast<int>(Nz_per_asmbly_) <= k) {
    return false;
  }

  // get the origin based on the fuelpin index
  const double origin_x0 =
      assembly_x0 + (static_cast<double>(i) + 0.5) * pitch_x_;
  const double origin_y0 =
      assembly_y0 + (static_cast<double>(j) + 0.5) * pitch_y_;

  const double xp = p.r.x() - origin_x0;
  const double yp = p.r.y() - origin_y0;
  const double radius_square = xp * xp + yp * yp;

  // Get the index of radial rings through linear search
  int r = -1;
  for (std::size_t ir = 0; ir < Nr_; ir++) {
    if (radius_square <= radius_ring_square_[ir]) {
      r = static_cast<int>(ir);
      break;
    }
  }

  // if the index of the rings is out of bounds, then don't keep the particle
  if (r < 0 || static_cast<int>(Nr_) <= r) {
    return false;
  }

  // Get the index of the angular segment
  int s = -1;
  if (xp >= 0.) {
    if (yp >= 0.) {
      s = 0;
    } else {
      s = 3;
    }
  } else {
    if (yp >= 0.) {
      s = 1;
    } else {
      s = 2;
    }
  }

  // get the key which will fit the particle into desired index
  std::uint32_t bin_key = static_cast<std::uint32_t>(
      l * Se_ + s * Ss_ + r * Sr_ + k * Sk_ + j * Sj_ + i * Si_);

  if (bins_.find(bin_key) == bins_.end()) {
    bins_[bin_key] = std::vector<BankedParticle*>();
  }

  bins_[bin_key].push_back(&p);

  return true;
}

std::vector<std::uint32_t> FuelPinApproxCancelator::sync_keys() {
  // Each node collects all keys
  std::vector<std::uint32_t> keys;
  keys.reserve(bins_.size());
  for (auto& key_bin_pair : bins_) {
    keys.push_back(key_bin_pair.first);
  }
  std::set<std::uint32_t> key_set;

  // Put master keys into the keyset
  if (mpi::rank == 0) {
    std::copy(keys.begin(), keys.end(), std::inserter(key_set, key_set.end()));
    keys.clear();
  }

  // For every Node starting at 1, send its keys to master and add to key_set
  for (int i = 1; i < mpi::size; i++) {
    if (mpi::rank == i) {
      auto nkeys = keys.size();
      mpi::Send(nkeys, 0);
      mpi::Send(std::span<std::uint32_t>(keys.begin(), keys.end()), 0);
      keys.clear();
    } else if (mpi::rank == 0) {
      std::size_t nkeys = 0;
      mpi::Recv(nkeys, i);
      keys.resize(nkeys);

      mpi::Recv(std::span<std::uint32_t>(keys.begin(), keys.end()), i);
      std::copy(keys.begin(), keys.end(),
                std::inserter(key_set, key_set.end()));
    }
  }

  // Move keys from key_set back to keys in master
  if (mpi::rank == 0) {
    keys.clear();
    keys.assign(key_set.begin(), key_set.end());
    key_set.clear();
  }

  // Send keys to all nodes from Master
  mpi::Bcast(keys, 0);

  return keys;
}

void FuelPinApproxCancelator::perform_cancellation_loop() {
  // Get keys of all non empty bins
  std::vector<std::uint32_t> keys = sync_keys();

  for (const auto key : keys) {
    std::uint64_t n_total = 0;
    double sum_wgt = 0.;
    double sum_wgt2 = 0.;

    // Go through all particles in the bin and add up positives and negatives
    // and get total sum
    auto& bin = bins_[key];
    for (const auto& p : bin) {
      n_total++;
      sum_wgt += p->wgt;
      sum_wgt2 += p->wgt2;
    }
    // Sum weights of all particles in each bin across all nodes
    mpi::Allreduce_sum(sum_wgt);
    if (this->cancel_dual_weights()) {
      mpi::Allreduce_sum(sum_wgt2);
    }

    // Sum total number of positive and negative particles across all nodes
    mpi::Allreduce_sum(n_total);

    // Set the avg weights
    const double inv_n = 1. / static_cast<double>(n_total);
    const double avg_wgt = sum_wgt * inv_n;
    const double avg_wgt2 = sum_wgt2 * inv_n;

    // Loop through particles in the bin once again and perform cancellation if
    // necessary
    for (const auto& p : bin) {
      p->wgt = avg_wgt;
      p->wgt2 = avg_wgt2;
    }
    bins_[key].clear();
  }
  keys.clear();
}

void FuelPinApproxCancelator::perform_cancellation_full_vector() {
  xt::xarray<double> wgts;
  const std::uint32_t Nx = Nx_asmbly_ * Nx_per_asmbly_;
  const std::uint32_t Ny = Ny_asmbly_ * Ny_per_asmbly_;
  const std::uint32_t Nz = Nz_asmbly_ * Nz_per_asmbly_;

  if (this->cancel_dual_weights()) {
    wgts.resize({2, Nx, Ny, Nz, Nr_, Ns_, Ne_});
  } else {
    wgts.resize({Nx, Ny, Nz, Nr_, Ns_, Ne_});
  }
  wgts.fill(0.);

  xt::xarray<std::uint32_t> n_totals;
  n_totals.resize({Nx, Ny, Nz, Nr_, Ns_, Ne_});
  n_totals.fill(0);

  for (auto& key_bin_pair : bins_) {
    std::uint32_t indx = key_bin_pair.first;
    const auto& bin = key_bin_pair.second;

    // Get i,j,k,r,s,l from the indx using strides
    std::uint32_t i = indx / Si_;
    std::uint32_t j = (indx - i * Si_) / Sj_;
    std::uint32_t k = (indx - i * Si_ - j * Sj_) / Sk_;
    std::uint32_t r = (indx - i * Si_ - j * Sj_ - k * Sk_) / Sr_;
    std::uint32_t s = (indx - i * Si_ - j * Sj_ - k * Sk_ - r * Sr_) / Ss_;
    std::uint32_t l =
        (indx - i * Si_ - j * Sj_ - k * Sk_ - r * Sr_ - s * Ss_) / Se_;

    std::uint32_t n_total = 0;
    double sum_wgt = 0.;
    double sum_wgt2 = 0.;

    // Go through all particles in the bin and add up positives and negatives
    // and get total sum
    for (const auto& p : bin) {
      n_total++;
      sum_wgt += p->wgt;
      sum_wgt2 += p->wgt2;
    }

    // Push the counts to the vectors
    n_totals(i, j, k, r, s, l) = n_total;

    if (this->cancel_dual_weights()) {
      wgts(0, i, j, k, r, s, l) = sum_wgt;
      wgts(1, i, j, k, r, s, l) = sum_wgt2;
    } else {
      wgts(i, j, k, r, s, l) = sum_wgt;
    }
  }

  // Sum the vectors across all nodes
  std::span<double> wgts_vals(wgts.data(), wgts.size());
  mpi::Allreduce_sum(wgts_vals);

  std::span<std::uint32_t> n_totals_vals(n_totals.data(), n_totals.size());
  mpi::Allreduce_sum(n_totals_vals);

  // all the vectors have size keys.size() so we use variable x to index them
  // since they should match to keys
  for (auto& key_bin_pair : bins_) {
    std::uint32_t indx = key_bin_pair.first;
    auto& bin = key_bin_pair.second;

    // Get i,j,k,r,s,l from the indx using strides
    std::uint32_t i = indx / Si_;
    std::uint32_t j = (indx - i * Si_) / Sj_;
    std::uint32_t k = (indx - i * Si_ - j * Sj_) / Sk_;
    std::uint32_t r = (indx - i * Si_ - j * Sj_ - k * Sk_) / Sr_;
    std::uint32_t s = (indx - i * Si_ - j * Sj_ - k * Sk_ - r * Sr_) / Ss_;
    std::uint32_t l =
        (indx - i * Si_ - j * Sj_ - k * Sk_ - r * Sr_ - s * Ss_) / Se_;

    // Set the avg weights
    const double inv_n = 1. / static_cast<double>(n_totals(i, j, k, r, s, l));
    double avg_wgt = 0.;
    double avg_wgt2 = 0.;
    if (this->cancel_dual_weights()) {
      avg_wgt = wgts(0, i, j, k, r, s, l) * inv_n;
      avg_wgt2 = wgts(1, i, j, k, r, s, l) * inv_n;
    } else {
      avg_wgt = wgts(i, j, k, r, s, l) * inv_n;
    }

    for (auto& p : bin) {
      p->wgt = avg_wgt;
      p->wgt2 = avg_wgt2;
    }

    bin.clear();
  }
}

void FuelPinApproxCancelator::perform_cancellation() {
  if (loop) {
    this->perform_cancellation_loop();
  } else {
    // this->perform_cancellation_vector();
    this->perform_cancellation_full_vector();
  }
}

std::shared_ptr<FuelPinApproxCancelator> make_fuelpin_approximate_cancelator(
    const YAML::Node& node) {
  // Get the low point
  if (!node["low"] || !node["low"].IsSequence() || (node["low"].size() != 3)) {
    fatal_error("No valid low entry for fuelpin approximate cancelator.");
  }
  double xl = node["low"][0].as<double>();
  double yl = node["low"][1].as<double>();
  double zl = node["low"][2].as<double>();
  Position r_low(xl, yl, zl);

  // Get the high point
  if (!node["high"] || !node["high"].IsSequence() ||
      (node["high"].size() != 3)) {
    fatal_error("No valid high entry for fuelpin approximate cancelator.");
  }
  double xh = node["high"][0].as<double>();
  double yh = node["high"][1].as<double>();
  double zh = node["high"][2].as<double>();
  Position r_high(xh, yh, zh);

  // Get the assemblies pitch in each direction
  if (!node["assembly-pitch"] || !node["assembly-pitch"].IsSequence() ||
      (node["assembly-pitch"].size() != 3)) {
    fatal_error(
        "No valid assembly-pitch entry for fuelpin approximate "
        "cancelator.");
  }
  std::vector<double> assembly_pitch =
      node["assembly-pitch"].as<std::vector<double>>();

  // Get the shape the assemblies in each direction
  if (!node["assembly-shape"] || !node["assembly-shape"].IsSequence() ||
      (node["assembly-shape"].size() != 3)) {
    fatal_error(
        "No valid assembly-shape entry for fuelpin approximate "
        "cancelator.");
  }
  std::vector<std::uint32_t> assembly_shape =
      node["assembly-shape"].as<std::vector<std::uint32_t>>();

  // get the fuelpin pitch
  if (!node["fuelpin-pitch"] || !node["fuelpin-pitch"].IsSequence() ||
      (node["fuelpin-pitch"].size() != 3)) {
    fatal_error(
        "No valid fuelpin-pitch entry for fuelpin approximate cancelator.");
  }
  std::vector<double> fuelpin_pitch =
      node["fuelpin-pitch"].as<std::vector<double>>();

  // Get the shape of the fuelpins per assemblies in each direction
  if (!node["fuelpin-per-assembly"] ||
      !node["fuelpin-per-assembly"].IsSequence() ||
      (node["fuelpin-per-assembly"].size() != 3)) {
    fatal_error(
        "No valid shape of fuelpin-per-assembly entry for fuelpin "
        "approximate cancelator.");
  }
  std::vector<std::uint32_t> fuelpins_per_assembly_shape =
      node["fuelpin-per-assembly"].as<std::vector<std::uint32_t>>();

  // get the inter-assembly gap
  if (!node["inter-assembly-gap"] || !node["inter-assembly-gap"].IsSequence() ||
      (node["inter-assembly-gap"].size() != 3)) {
    fatal_error(
        "No valid inter-assembly-gap entry for fuelpin "
        "approximate cancelator.");
  }
  std::vector<double> inter_assembly_gap =
      node["inter-assembly-gap"].as<std::vector<double>>();

  // get the axial orientation of the fuel-pin
  if (!node["axial-orientation"] || !node["axial-orientation"].IsScalar()) {
    fatal_error(
        "No valid axial-orientation entry for fuelpin approximate cancelator.");
  }

  std::string string_axis_orient = node["axial-orientation"].as<std::string>();
  FuelPinApproxCancelator::Orientation axis_orient;
  if (string_axis_orient == "x") {
    axis_orient = FuelPinApproxCancelator::Orientation::X;
  } else if (string_axis_orient == "y") {
    axis_orient = FuelPinApproxCancelator::Orientation::Y;
  } else if (string_axis_orient == "z") {
    axis_orient = FuelPinApproxCancelator::Orientation::Z;
  } else {
    fatal_error(
        "axial-orient can only be either x, y, or z for fuelpin approxiate "
        "cancellator.");
  }

  // check the entry for the loop
  bool loop = false;
  if (node["loop"] && node["loop"].IsScalar()) {
    loop = node["loop"].as<bool>();
  } else if (node["loop"]) {
    fatal_error("Invalid entry for loop in fuelpin approximate cancelator.");
  }

  // check the entry for the energy-bounds.
  std::vector<double> energy_bounds;
  if (node["energy-bounds"] && node["energy-bounds"].IsSequence()) {
    energy_bounds = node["energy-bounds"].as<std::vector<double>>();
  } else if (node["energy-bounds"]) {
    fatal_error(
        "No valid energy-bounds entry for fuelpin approximate cancelator.");
  }

  // check which method is used to divide the disc into different rings.
  bool is_N_ring_given = false;
  if (node["radius"] && node["rings"]) {
    is_N_ring_given = true;
  }

  bool is_radii_given = false;
  if (node["radii"]) {
    is_radii_given = true;
  }

  if (is_N_ring_given && is_radii_given) {
    fatal_error(
        "Two methods, equi-volume-rings and rings radii are given, cannot be "
        "used at the same time in fuelpin approximate cancelator.");
  }
  if ((!is_N_ring_given) && (!is_radii_given)) {
    fatal_error(
        "Neither equi-volume-rings nor rings radii are given for fuelpin "
        "approximate cancelator.");
  }

  double radius;
  std::uint32_t Nr;
  if (is_N_ring_given) {
    if (!node["radius"].IsScalar()) {
      fatal_error(
          "No valid radius entry is given for the fuelpin approximate "
          "cancelator.");
    }
    if (!node["rings"].IsScalar()) {
      fatal_error(
          "No valid rings entry is given for the fuelpin approximate "
          "cancelator.");
    }
    radius = node["radius"].as<double>();
    Nr = node["rings"].as<std::uint32_t>();
  }

  std::vector<double> radii;
  if (is_radii_given) {
    if (!node["radii"].IsSequence()) {
      fatal_error("No valid radii entry is given for the fuelpin approximate.");
    }
    radii = node["radii"].as<std::vector<double>>();
  }

  if (loop) {
    Output::instance().write(" Using FuelPinApproxCancelator with loop.\n");
  } else {
    Output::instance().write(" Using FuelPinApproxCancelator with vector.\n");
  }

  // create the FuelPinApproxCancelator based on the
  if (is_N_ring_given) {
    // construct based on the number of rings and radius of the cylinder
    if (energy_bounds.size() == 0) {
      return std::make_shared<FuelPinApproxCancelator>(
          r_low, r_high, assembly_pitch, assembly_shape, fuelpin_pitch,
          fuelpins_per_assembly_shape, inter_assembly_gap, radius, Nr,
          axis_orient, loop);
    } else {
      return std::make_shared<FuelPinApproxCancelator>(
          r_low, r_high, assembly_pitch, assembly_shape, fuelpin_pitch,
          fuelpins_per_assembly_shape, inter_assembly_gap, radius, Nr,
          axis_orient, energy_bounds, loop);
    }
  }

  // otherwise construct based on the radii of the different rings
  if (energy_bounds.size() == 0) {
    return std::make_shared<FuelPinApproxCancelator>(
        r_low, r_high, assembly_pitch, assembly_shape, fuelpin_pitch,
        fuelpins_per_assembly_shape, inter_assembly_gap, radii, axis_orient,
        loop);
  } else {
    return std::make_shared<FuelPinApproxCancelator>(
        r_low, r_high, assembly_pitch, assembly_shape, fuelpin_pitch,
        fuelpins_per_assembly_shape, inter_assembly_gap, radii, axis_orient,
        energy_bounds, loop);
  }
}