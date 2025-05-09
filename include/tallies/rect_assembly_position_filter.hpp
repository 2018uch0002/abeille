#ifndef RECT_ASSEMBLY_POSITION_FILTER_H
#define RECT_ASSEMBLY_POSITION_FILTER_H

#include <tallies/cartesian_filter.hpp>
#include <utils/position.hpp>

#include <array>

// RectAssemblyPositionFilter
class RectAssemblyPositionFilter : public CartesianFilter {
 public:
  RectAssemblyPositionFilter(Position r_low, Position r_high, std::vector<std::size_t> assembly_shape,
                             std::vector<double> inter_assembly_gap, std::vector<std::size_t> bin_shape_per_assembly, 
                             std::size_t id);

  StaticVector3 get_indices(const Tracker& tktr) const override final;

  StaticVector3 get_position_index(const Position& r) const override final;

  std::vector<TracklengthDistance> get_indices_tracklength(
      const Tracker& /*trkr*/, double /*d_flight*/) const override final{
    
    fatal_error("track-length capabilites for the react-assembly-position-filter is not implemented yet.");

    std::vector<TracklengthDistance> trk_len_distance;
    return trk_len_distance;
  }

  double x_min(const StaticVector3& index) const override final;
  double x_max(const StaticVector3& index) const override final;
  double dx(const StaticVector3& /*index*/) const override final { return dx_bin_; }
  double inv_dx(const StaticVector3& /*index*/) const override final {
    return inv_dx_bin_;
  }
  double y_min(const StaticVector3& index) const override final;
  double y_max(const StaticVector3& index) const override final;
  double dy(const StaticVector3& /*index*/) const override final { return dy_bin_; }
  double inv_dy(const StaticVector3& /*index*/) const override final {
    return inv_dy_bin_;
  }
  double z_min(const StaticVector3& index) const override final;
  double z_max(const StaticVector3& index) const override final;
  double dz(const StaticVector3& /*index*/) const override final { return dz_bin_; }
  double inv_dz(const StaticVector3& /*index*/) const override final {
    return inv_dz_bin_;
  }

  StaticVector3 get_shape() const override final {
    if (real_Nx_ == 1 && real_Ny_ == 1 && real_Nz_ == 1) {
      return {1};
    }
    return reduce_dimension(real_Nx_, real_Ny_, real_Nz_);
  }

  std::string type_str() const override { return "rect-assembly-position-filter"; }

 protected:
  // required for track-length
  // bool find_entry_point(Position& r, const Direction& u,
  //                       double& d_flight) const;
  // void initialize_indices(const Position& r, const Direction& u, int& i, int& j,
  //                         int& k, std::array<int, 3>& on) const;
  // void update_indices(int key, int& i, int& j, int& k,
  //                     std::array<int, 3>& on) const;

  // std::pair<double, int> distance_to_next_index(const Position& r,
  //                                               const Direction& u,
  //                                               const std::array<int, 3>& on,
  //                                               int i, int j, int k) const;

  void write_to_hdf5(H5::Group& grp) const override final;

 private:
  double assembly_dx_, assembly_dy_, assembly_dz_; // assembly pitches
  double inv_assembly_dx_, inv_assembly_dy_, inv_assembly_dz_; // inverse of assembly pitches
  double inter_asmbly_gap_x_, inter_asmbly_gap_y_, inter_asmbly_gap_z_; // inter-assembly gaps 
  double dx_bin_, dy_bin_, dz_bin_; // width of a bin in x, y, and z
  double inv_dx_bin_, inv_dy_bin_, inv_dz_bin_; // inverse of width of a bin in x, y, and z  
  
  std::size_t asmbly_Nx_, asmbly_Ny_, asmbly_Nz_; // shape of the assembly
  std::size_t Nx_, Ny_, Nz_; // shape of the bin per assembly
  std::size_t N_gap_x_, N_gap_y_, N_gap_z_; // distrization of gaps in x, y, and z on one side 
  std::size_t real_Nx_, real_Ny_, real_Nz_, x_index_, y_index_, z_index_; // real number of the bins in the filter

  // function will reduce the dimsion, if there is only one bin in the direction
  StaticVector3 reduce_dimension(const size_t& loc_x, const size_t& loc_y,
                                 const size_t& loc_z) const {
    if (real_Nx_ == 1 && real_Ny_ == 1 && real_Nz_ == 1) {
      return {loc_x};
    }
    StaticVector3 reduce_;
    if (real_Nx_ > 1) {
      reduce_.push_back(loc_x);
    }

    if (real_Ny_ > 1) {
      reduce_.push_back(loc_y);
    }

    if (real_Nz_ > 1) {
      reduce_.push_back(loc_z);
    }
    return reduce_;
  }
};

// Make the cartesian or position filter class
std::shared_ptr<RectAssemblyPositionFilter> make_rect_assembly_position_filter(
    const YAML::Node& node);

#endif