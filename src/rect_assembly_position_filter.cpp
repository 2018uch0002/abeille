#include <tallies/rect_assembly_position_filter.hpp>
#include <utils/constants.hpp>
#include <utils/output.hpp>

#include <sstream>

RectAssemblyPositionFilter::RectAssemblyPositionFilter(Position r_low, Position r_high, 
      std::vector<std::size_t> assembly_shape,
      std::vector<double> inter_assembly_gap, std::vector<std::size_t> bin_shape_per_assembly, 
      std::size_t id)
    : CartesianFilter(r_low, r_high, id),
      assembly_dx_(),
      assembly_dy_(),
      assembly_dz_(),
      inv_assembly_dx_(),
      inv_assembly_dy_(),
      inv_assembly_dz_(),
      inter_asmbly_gap_x_(),
      inter_asmbly_gap_y_(),
      inter_asmbly_gap_z_(),
      dx_bin_(),
      dy_bin_(),
      dz_bin_(),
      inv_dx_bin_(),
      inv_dy_bin_(),
      inv_dz_bin_(),
      asmbly_Nx_(),
      asmbly_Ny_(),
      asmbly_Nz_(),
      Nx_(),
      Ny_(),
      Nz_(),
      real_Nx_(),
      real_Ny_(),
      real_Nz_(),
      x_index_(),
      y_index_(),
      z_index_() {
  if ((r_low_.x() >= r_high_.x()) || (r_low_.y() >= r_high_.y()) ||
      (r_low_.z() >= r_high_.z())){
    fatal_error("In position-filter with id: " + std::to_string(id) +
                ", coordinates of \"low\" position are >= than \"high\" "
                "position.");
  }

  // check the shape of the assembly
  if (assembly_shape.size() != 3){
    fatal_error("In position-filter with id: " + std::to_string(id) +
                ", the shape of the assembly must be of size 3.");
  }

  asmbly_Nx_ = assembly_shape[0];
  asmbly_Ny_ = assembly_shape[1];
  asmbly_Ny_ = assembly_shape[2];
  if (asmbly_Nx_ == 0 || asmbly_Ny_ == 0 || asmbly_Nz_ == 0){
    fatal_error("In position-filter with id: " + std::to_string(id) + ", the number of bins in any direction must be a non-zero.");
  }

  assembly_dx_ = (r_high_.x() - r_low_.x()) / static_cast<double>(asmbly_Nx_);
  assembly_dy_ = (r_high_.y() - r_low_.y()) / static_cast<double>(asmbly_Ny_);
  assembly_dz_ = (r_high_.z() - r_low_.z()) / static_cast<double>(asmbly_Nz_);

  inv_assembly_dx_ = 1./ assembly_dx_;
  inv_assembly_dy_ = 1./ assembly_dy_;
  inv_assembly_dz_ = 1./ assembly_dz_;

  // check the shape of the bins per assembly
  if (bin_shape_per_assembly.size() != 3){
    fatal_error("In position-filter with id: " + std::to_string(id) + ", the shape of the bin per assembly must be of size 3.");
  }

  Nx_ = bin_shape_per_assembly[0];
  Ny_ = bin_shape_per_assembly[1];
  Nz_ = bin_shape_per_assembly[2];

  if (Nx_ == 0 || Ny_ == 0 || Nz_ == 0){
    fatal_error("In position-filter with id: " + std::to_string(id) + ", the number of bins per assembly in any direction must be a non-zero.");
  }

  // get the inter assembly_gaps
  if (inter_assembly_gap.size() != 3){
    fatal_error("In position-filter with id: " + std::to_string(id) + ", the size of inter-assembly-gap must be 3.");
  }

  inter_asmbly_gap_x_ = inter_assembly_gap[0];
  inter_asmbly_gap_y_ = inter_assembly_gap[1];
  inter_asmbly_gap_z_ = inter_assembly_gap[2];

  if (inter_asmbly_gap_x_ < 0. || inter_asmbly_gap_y_ < 0. || inter_asmbly_gap_z_ < 0.){
    fatal_error("In position-filter with id: " + std::to_string(id) + ", the inter-assembly-gap should not be a negative number.");
  }

  dx_bin_ = (assembly_dx_ - inter_asmbly_gap_x_) / static_cast<double>(Nx_);
  dy_bin_ = (assembly_dy_ - inter_asmbly_gap_y_) / static_cast<double>(Ny_);
  dz_bin_ = (assembly_dz_ - inter_asmbly_gap_z_) / static_cast<double>(Nz_); 

  inv_dx_bin_ = 1. / dx_bin_;
  inv_dy_bin_ = 1. / dy_bin_;
  inv_dz_bin_ = 1. / dz_bin_;

  real_Nx_ = Nx_ * asmbly_Nx_;
  real_Ny_ = Ny_ * asmbly_Ny_;
  real_Nz_ = Nz_ * asmbly_Nz_;

  x_index_ = 0;
  y_index_ = 1;
  z_index_ = 2;
  if (real_Nx_ == 1) {
    x_index_ = 0;
    y_index_--;
    z_index_--;
  }

  if (real_Ny_ == 1) {
    y_index_ = 0;
    z_index_--;
  }

  if (real_Nz_ == 1) {
    z_index_ = 0;
  }
}

StaticVector3 RectAssemblyPositionFilter::get_indices(
    const Tracker& tktr) const {
  StaticVector3 indices;
  const Position r = tktr.r();

  const int index_asmbly_x = static_cast<int>(std::floor((r.x() - r_low_.x()) * inv_assembly_dx_));
  const int index_asmbly_y = static_cast<int>(std::floor((r.y() - r_low_.y()) * inv_assembly_dy_));
  const int index_asmbly_z = static_cast<int>(std::floor((r.z() - r_low_.z()) * inv_assembly_dz_));
  
  if ((index_asmbly_x >=0  && index_asmbly_x < static_cast<int>(asmbly_Nx_)) &&
      (index_asmbly_y >=0  && index_asmbly_y < static_cast<int>(asmbly_Ny_)) &&
      (index_asmbly_z >=0  && index_asmbly_z < static_cast<int>(asmbly_Nz_))){
    
        // if we are here, means we are in some assembly

        const double bin_xmin = r_low_.x() + assembly_dx_ * index_asmbly_x + inter_asmbly_gap_x_ * 0.5;
        const double bin_ymin = r_low_.y() + assembly_dy_ * index_asmbly_y + inter_asmbly_gap_y_ * 0.5;
        const double bin_zmin = r_low_.z() + assembly_dz_ * index_asmbly_z + inter_asmbly_gap_z_ * 0.5;
        
        const int index_bin_x = static_cast<int>(std::floor((r.x() - bin_xmin) * inv_dx_bin_));
        const int index_bin_y = static_cast<int>(std::floor((r.y() - bin_ymin) * inv_dy_bin_));
        const int index_bin_z = static_cast<int>(std::floor((r.z() - bin_zmin) * inv_dz_bin_));

        if ((index_bin_x >= 0 && index_bin_x < static_cast<int>(Nx_)) &&
            (index_bin_y >= 0 && index_bin_y < static_cast<int>(Ny_)) &&
            (index_bin_z >= 0 && index_bin_z < static_cast<int>(Nz_))){
              indices = reduce_dimension(index_asmbly_x * Nx_ + index_bin_x, index_asmbly_y * Ny_ + index_bin_y, index_asmbly_z * Nz_ + index_bin_z);
            }
  }

  return indices;
}

StaticVector3 RectAssemblyPositionFilter::get_position_index(
    const Position& r) const {
      fatal_error("yet to be implemented");
  // StaticVector3 indices;
  // const int index_x =
  //     static_cast<int>(std::floor((r.x() - r_low_.x()) * dx_inv_));
  // const int index_y =
  //     static_cast<int>(std::floor((r.y() - r_low_.y()) * dy_inv_));
  // const int index_z =
  //     static_cast<int>(std::floor((r.z() - r_low_.z()) * dz_inv_));

  // if ((index_x >= 0 && index_x < static_cast<int>(Nx_)) &&
  //     (index_y >= 0 && index_y < static_cast<int>(Ny_)) &&
  //     (index_z >= 0 && index_z < static_cast<int>(Nz_))) {
  //   indices = reduce_dimension(index_x, index_y, index_z);
  // }

  // return indices;
}

double RectAssemblyPositionFilter::x_min(const StaticVector3& index) const {
  if (real_Nx_ == 1) {
    return r_low_.x();
  }
  const int index_assmbly_x = static_cast<int>(std::floor(index[x_index_] / Nx_));
  const int index_bin_x = index[x_index_] - index_assmbly_x * Nx_;
  const double xmin_bin = r_low_.x() + static_cast<double>(index_assmbly_x) * assembly_dx_ + 0.5 * inter_asmbly_gap_x_ + static_cast<double>(index_bin_x) * dx_bin_;
  return xmin_bin;
}

double RectAssemblyPositionFilter::x_max(const StaticVector3& index) const {
  if (real_Nx_ == 1) {
    return r_high_.x();
  }
  const int index_assmbly_x = static_cast<int>(std::floor(index[x_index_] / Nx_));
  const int index_bin_x = index[x_index_] - index_assmbly_x * Nx_;
  const double xmax_bin = r_low_.x() + static_cast<double>(index_assmbly_x) * assembly_dx_ + 0.5 * inter_asmbly_gap_x_ + (static_cast<double>(index_bin_x) + 1) * dx_bin_;
  
  return xmax_bin;
}

double RectAssemblyPositionFilter::y_min(const StaticVector3& index) const {
  if (real_Ny_ == 1) {
    return r_low_.y();
  }

  const int index_assmbly_y = static_cast<int>(std::floor(index[y_index_] / Ny_));
  const int index_bin_y = index[y_index_] - index_assmbly_y * Ny_;
  const double ymin_bin = r_low_.y() + static_cast<double>(index_assmbly_y) * assembly_dy_ + 0.5 * inter_asmbly_gap_y_ + static_cast<double>(index_bin_y) * dy_bin_;
  
  return ymin_bin;
}

double RectAssemblyPositionFilter::y_max(const StaticVector3& index) const {
  if (real_Ny_ == 1) return r_high_.y();
  
  const int index_assmbly_y = static_cast<int>(std::floor(index[y_index_] / Ny_));
  const int index_bin_y = index[y_index_] - index_assmbly_y * Ny_;
  const double ymax_bin = r_low_.y() + static_cast<double>(index_assmbly_y) * assembly_dy_ + 0.5 * inter_asmbly_gap_y_ + (static_cast<double>(index_bin_y) + 1) * dy_bin_;
  
  return ymax_bin;
}

double RectAssemblyPositionFilter::z_min(const StaticVector3& index) const {
  if (real_Nz_ == 1) return r_low_.z();
 
  const int index_assmbly_z = static_cast<int>(std::floor(index[z_index_] / Nz_));
  const int index_bin_z = index[z_index_] - index_assmbly_z * Nz_;
  const double zmin_bin = r_low_.z() + static_cast<double>(index_assmbly_z) * assembly_dz_ + 0.5 * inter_asmbly_gap_z_ + static_cast<double>(index_bin_z) * dz_bin_;

  return zmin_bin;
}

double RectAssemblyPositionFilter::z_max(const StaticVector3& index) const {
  if (real_Nz_ == 1) return r_high_.z();

  const int index_assmbly_z = static_cast<int>(std::floor(index[z_index_] / Nz_));
  const int index_bin_z = index[z_index_] - index_assmbly_z * Nz_;
  const double zmax_bin = r_low_.z() + static_cast<double>(index_assmbly_z) * assembly_dz_ + 0.5 * inter_asmbly_gap_z_ + (static_cast<double>(index_bin_z) + 1)* dz_bin_;

  return zmax_bin;
}




void RectAssemblyPositionFilter::write_to_hdf5(H5::Group& grp) const {
  // Save id in attributes
  if (grp.hasAttribute("id")) {
    grp.deleteAttribute("id");
  }
  grp.createAttribute("id", this->id());

  // Save type in attributes
  if (grp.hasAttribute("type")) {
    grp.deleteAttribute("type");
  }
  grp.createAttribute("type", "regular-cartesian-mesh");

  // Save low position
  std::array<double, 3> r_low{r_low_.x(), r_low_.y(), r_low_.z()};
  if (grp.hasAttribute("low")) {
    grp.deleteAttribute("low");
  }
  grp.createAttribute("low", r_low);

  // Save high position
  std::array<double, 3> r_high{r_high_.x(), r_high_.y(), r_high_.z()};
  if (grp.hasAttribute("high")) {
    grp.deleteAttribute("high");
  }
  grp.createAttribute("high", r_high);

  // Save shape
  std::array<std::size_t, 3> shape{Nx_, Ny_, Nz_};
  if (grp.hasAttribute("shape")) {
    grp.deleteAttribute("shape");
  }
  grp.createAttribute("shape", shape);

  std::vector<double> x_bounds(real_Nx_ + 1, 0.);
  std::size_t itr = 0;
  for (std::size_t i = 0; i < asmbly_Nx_; i++ ){
    double x0 = r_low_.x() + i * assembly_dx_ + inter_asmbly_gap_x_ * 0.5;
    for (std::size_t j = 0; j <= Nx_; j++){
      x_bounds[itr] = x0 + dx_bin_ * j;
      itr++;
    }
  }
  grp.createDataSet("x-bounds", x_bounds);

  std::vector<double> y_bounds(real_Ny_ + 1, 0.);
  itr = 0;
  for (std::size_t i = 0; i < asmbly_Ny_; i++ ){
    double y0 = r_low_.y() + i * assembly_dy_ + inter_asmbly_gap_y_ * 0.5;
    for (std::size_t j = 0; j <= Ny_; j++){
      y_bounds[itr] = y0 + dy_bin_ * j;
      itr++;
    }
  }
  grp.createDataSet("y-bounds", y_bounds);

  std::vector<double> z_bounds(real_Nz_ + 1, 0.);
  itr = 0;
  for (std::size_t i = 0; i < asmbly_Nz_; i++ ){
    double z0 = r_low_.z() + i * assembly_dz_ + inter_asmbly_gap_z_ * 0.5;
    for (std::size_t j = 0; j <= Nz_; j++){
      z_bounds[itr] = z0 + dz_bin_ * j;
      itr++;
    }
  }
  grp.createDataSet("z-bounds", z_bounds);
}

// Make the cartesian or position filter class
std::shared_ptr<RectAssemblyPositionFilter> make_rect_assembly_position_filter(
    const YAML::Node& node) {
  if (!node["id"] || !node["id"].IsScalar()) {
    fatal_error("Invalid id is given for the position-filter.");
  }
  std::size_t id = node["id"].as<std::size_t>();

  if (!node["low"]) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", \"low\" coordinates are not provided.";
    fatal_error(mssg.str());
  } else if (!node["low"].IsSequence() || node["low"].size() != 3) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", the given entry for the \"low\" coordinates must be a sequence "
            "of size 3.";
    fatal_error(mssg.str());
  }

  if (!node["high"]) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", \"high\" coordinates are not provided.";
    fatal_error(mssg.str());
  } else if (!node["high"].IsSequence() || node["high"].size() != 3) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", the given entry for the \"high\" coordinates must be a sequence "
            "of size 3.";
    fatal_error(mssg.str());
  }

  // check for the shape of the assemblies
  if (!node["assembly-shape"]) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
          << ", \"assembly-shape\" is not provided.";
    fatal_error(mssg.str());
  } else if (!node["assembly-shape"].IsSequence() || node["assembly-shape"].size() != 3) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
          << ", \"assembly-shape\" must be a sequence of size 3.";
    fatal_error(mssg.str());
  }

  // check for the shape per assembly
  if (!node["bin-shape-per-assembly"]) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", \"bin-shape-per-assembly\" is not provided.";
    fatal_error(mssg.str());
  } else if (!node["bin-shape-per-assembly"].IsSequence() || node["bin-shape-per-assembly"].size() != 3) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", \"bin-shape-per-assembly\" must be a sequence of size 3.";
    fatal_error(mssg.str());
  }

  // get the inter-assembly-gap
  if (!node["inter-assembly-gap"]){
    std::stringstream mssg;
    mssg << "For position-filter with id " << id << ", \"inter-assembly-gaps\" are not provided.";
    fatal_error(mssg.str());
  } else if (!!node["inter-assembly-gap"].IsSequence() || node["inter-assembly-gap"].size() != 3){
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", \"inter-assembly-gap\" must be a sequence of size 3.";
    fatal_error(mssg.str());
  }

  std::vector<double> low_point = node["low"].as<std::vector<double>>();
  std::vector<double> high_point = node["high"].as<std::vector<double>>();
  std::vector<std::size_t> assembly_shape = node["assembly-shape"].as<std::vector<std::size_t>>();
  std::vector<std::size_t> bin_shape_per_asmbly = node["bin-shape-per-assembly"].as<std::vector<std::size_t>>();
  std::vector<double> inter_assembly_gap = node["inter-assembly-gap"].as<std::vector<double>>();

  Position r_low(low_point[0], low_point[1], low_point[2]);
  Position r_high(high_point[0], high_point[1], high_point[2]);

  std::shared_ptr<RectAssemblyPositionFilter> mesh_type_filter =
      std::make_shared<RectAssemblyPositionFilter>(r_low, r_high, assembly_shape, inter_assembly_gap, bin_shape_per_asmbly, id);

  return mesh_type_filter;
}
