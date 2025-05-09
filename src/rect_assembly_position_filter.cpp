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
      N_gap_x_(),
      N_gap_y_(),
      N_gap_z_(),
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
  asmbly_Nz_ = assembly_shape[2];
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

  if (inter_asmbly_gap_x_ == 0.) N_gap_x_ = 0;
  else N_gap_x_ = 1;

  if (inter_asmbly_gap_y_ == 0.) N_gap_y_ = 0;
  else N_gap_y_ = 1;

  if (inter_asmbly_gap_z_ == 0.) N_gap_z_ = 0;
  else N_gap_z_ = 1;

  dx_bin_ = (assembly_dx_ - inter_asmbly_gap_x_) / static_cast<double>(Nx_);
  dy_bin_ = (assembly_dy_ - inter_asmbly_gap_y_) / static_cast<double>(Ny_);
  dz_bin_ = (assembly_dz_ - inter_asmbly_gap_z_) / static_cast<double>(Nz_); 

  inv_dx_bin_ = 1. / dx_bin_;
  inv_dy_bin_ = 1. / dy_bin_;
  inv_dz_bin_ = 1. / dz_bin_;

  // since the gaps are usually very small therefore, 
  // the gaps will not be divided into multiple bins on the normal directions
  real_Nx_ = (N_gap_x_ + Nx_ + N_gap_x_) * asmbly_Nx_;
  real_Ny_ = (N_gap_y_ + Ny_ + N_gap_y_) * asmbly_Ny_;
  real_Nz_ = (N_gap_z_ + Nz_ + N_gap_z_) * asmbly_Nz_;

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

      const Position r = tktr.r();
      StaticVector3 indices;

      const int index_asmbly_x = static_cast<int>(std::floor((r.x() - r_low_.x()) * inv_assembly_dx_));
      const int index_asmbly_y = static_cast<int>(std::floor((r.y() - r_low_.y()) * inv_assembly_dy_));
      const int index_asmbly_z = static_cast<int>(std::floor((r.z() - r_low_.z()) * inv_assembly_dz_));
      
      if ((index_asmbly_x >=0  && index_asmbly_x < static_cast<int>(asmbly_Nx_)) &&
          (index_asmbly_y >=0  && index_asmbly_y < static_cast<int>(asmbly_Ny_)) &&
          (index_asmbly_z >=0  && index_asmbly_z < static_cast<int>(asmbly_Nz_))){
        
            // if we are here, means we are some here inside the assembly.
    
            const double bin_xmin = r_low_.x() + assembly_dx_ * index_asmbly_x + inter_asmbly_gap_x_ * 0.5;
            const double bin_ymin = r_low_.y() + assembly_dy_ * index_asmbly_y + inter_asmbly_gap_y_ * 0.5;
            const double bin_zmin = r_low_.z() + assembly_dz_ * index_asmbly_z + inter_asmbly_gap_z_ * 0.5;
            
            const int index_bin_x = static_cast<int>(std::floor((r.x() - bin_xmin) * inv_dx_bin_));
            const int index_bin_y = static_cast<int>(std::floor((r.y() - bin_ymin) * inv_dy_bin_));
            const int index_bin_z = static_cast<int>(std::floor((r.z() - bin_zmin) * inv_dz_bin_));
           
             
            std::size_t ix, iy, iz;
            
            // check if you are in the region of the excluding the gaps along x-direction 
            if (index_bin_x >= 0 && index_bin_x < static_cast<int>(Nx_)){
              ix = index_asmbly_x * (N_gap_x_+Nx_+N_gap_x_) + (N_gap_x_ + index_bin_x);
            } else if (index_bin_x < 0){ // you are in the gap regions along x-directions
              ix = index_asmbly_x * (N_gap_x_+Nx_+N_gap_x_) + (N_gap_x_-1);
            } else if (index_bin_x >= static_cast<int>(Nx_)){
              ix = index_asmbly_x * (N_gap_x_+Nx_+N_gap_x_) + (N_gap_x_ + Nx_+ N_gap_x_ - 1);
            } else {
              fatal_error("postion-filter id: " + std::to_string(this->id()) + ", it should not come here.");
            }
            
            // check if you are in the region of the excluding the gaps along y-direction 
            if (index_bin_y >= 0 && index_bin_y < static_cast<int>(Ny_)){
              iy = index_asmbly_y * (N_gap_y_+Ny_+N_gap_y_) + (N_gap_y_ + index_bin_y);
            } else if (index_bin_y < 0){ // you are in the gap regions along y-directions
              iy = index_asmbly_y * (N_gap_y_+Ny_+N_gap_y_) + (N_gap_y_-1);
            } else if (index_bin_y >= static_cast<int>(Ny_)){
              iy = index_asmbly_y * (N_gap_y_+Ny_+N_gap_y_) + (N_gap_y_ + Ny_+ N_gap_y_ - 1);
            } else {
              fatal_error("postion-filter id: " + std::to_string(this->id()) + ", it should not come here.");
            }
    
            // check if you are in the region of the excluding the gaps along z-direction 
            if (index_bin_z >= 0 && index_bin_z < static_cast<int>(Nz_)){
              iz = index_asmbly_z * (N_gap_z_+Nz_+N_gap_z_) + N_gap_z_ + index_bin_z;
            } else if (index_bin_z < 0){ // you are in the gap regions along x-directions
              ix = index_asmbly_z * (N_gap_z_+Nz_+N_gap_z_) + (N_gap_z_-1);
            } else if (index_bin_z >= static_cast<int>(Nz_)){
              ix = index_asmbly_z * (N_gap_z_+Nz_+N_gap_z_) + (N_gap_z_ + Nz_ + N_gap_z_ - 1);
            } else {
              fatal_error("postion-filter id: " + std::to_string(this->id()) + ", it should not come here.");
            }
    
            indices = reduce_dimension(ix, iy, iz);
      }

  return indices;
}

StaticVector3 RectAssemblyPositionFilter::get_position_index(
    const Position& r) const {

  StaticVector3 indices;

  const int index_asmbly_x = static_cast<int>(std::floor((r.x() - r_low_.x()) * inv_assembly_dx_));
  const int index_asmbly_y = static_cast<int>(std::floor((r.y() - r_low_.y()) * inv_assembly_dy_));
  const int index_asmbly_z = static_cast<int>(std::floor((r.z() - r_low_.z()) * inv_assembly_dz_));
  
  if ((index_asmbly_x >=0  && index_asmbly_x < static_cast<int>(asmbly_Nx_)) &&
      (index_asmbly_y >=0  && index_asmbly_y < static_cast<int>(asmbly_Ny_)) &&
      (index_asmbly_z >=0  && index_asmbly_z < static_cast<int>(asmbly_Nz_))){
    
        // if we are here, means we are some here inside the assembly.

        const double bin_xmin = r_low_.x() + assembly_dx_ * index_asmbly_x + inter_asmbly_gap_x_ * 0.5;
        const double bin_ymin = r_low_.y() + assembly_dy_ * index_asmbly_y + inter_asmbly_gap_y_ * 0.5;
        const double bin_zmin = r_low_.z() + assembly_dz_ * index_asmbly_z + inter_asmbly_gap_z_ * 0.5;
        
        const int index_bin_x = static_cast<int>(std::floor((r.x() - bin_xmin) * inv_dx_bin_));
        const int index_bin_y = static_cast<int>(std::floor((r.y() - bin_ymin) * inv_dy_bin_));
        const int index_bin_z = static_cast<int>(std::floor((r.z() - bin_zmin) * inv_dz_bin_));
       
         
        std::size_t ix, iy, iz;
        
        // check if you are in the region of the excluding the gaps along x-direction 
        if (index_bin_x >= 0 && index_bin_x < static_cast<int>(Nx_)){
          ix = index_asmbly_x * (N_gap_x_+Nx_+N_gap_x_) + (N_gap_x_ + index_bin_x);
        } else if (index_bin_x < 0){ // you are in the gap regions along x-directions
          ix = index_asmbly_x * (N_gap_x_+Nx_+N_gap_x_) + (N_gap_x_-1);
        } else if (index_bin_x >= static_cast<int>(Nx_)){
          ix = index_asmbly_x * (N_gap_x_+Nx_+N_gap_x_) + (N_gap_x_ + Nx_+ N_gap_x_ - 1);
        } else {
          fatal_error("postion-filter id: " + std::to_string(this->id()) + ", it should not come here.");
        }
        
        // check if you are in the region of the excluding the gaps along y-direction 
        if (index_bin_y >= 0 && index_bin_y < static_cast<int>(Ny_)){
          iy = index_asmbly_y * (N_gap_y_+Ny_+N_gap_y_) + (N_gap_y_ + index_bin_y);
        } else if (index_bin_y < 0){ // you are in the gap regions along y-directions
          iy = index_asmbly_y * (N_gap_y_+Ny_+N_gap_y_) + (N_gap_y_-1);
        } else if (index_bin_y >= static_cast<int>(Ny_)){
          iy = index_asmbly_y * (N_gap_y_+Ny_+N_gap_y_) + (N_gap_y_ + Ny_+ N_gap_y_ - 1);
        } else {
          fatal_error("postion-filter id: " + std::to_string(this->id()) + ", it should not come here.");
        }

        // check if you are in the region of the excluding the gaps along z-direction 
        if (index_bin_z >= 0 && index_bin_z < static_cast<int>(Nz_)){
          iz = index_asmbly_z * (N_gap_z_+Nz_+N_gap_z_) + N_gap_z_ + index_bin_z;
        } else if (index_bin_z < 0){ // you are in the gap regions along x-directions
          ix = index_asmbly_z * (N_gap_z_+Nz_+N_gap_z_) + (N_gap_z_-1);
        } else if (index_bin_z >= static_cast<int>(Nz_)){
          ix = index_asmbly_z * (N_gap_z_+Nz_+N_gap_z_) + (N_gap_z_ + Nz_ + N_gap_z_ - 1);
        } else {
          fatal_error("postion-filter id: " + std::to_string(this->id()) + ", it should not come here.");
        }

        indices = reduce_dimension(ix, iy, iz);
  }
  return indices;
}

double RectAssemblyPositionFilter::x_min(const StaticVector3& index) const {
  if (real_Nx_ == 1) {
    return r_low_.x();
  }

  const std::size_t index_assmbly_x = static_cast<std::size_t>(std::floor(index[x_index_] / (N_gap_x_ + Nx_ + N_gap_x_)));
  const std::size_t index_bin_x = static_cast<std::size_t>(index[x_index_] - index_assmbly_x * (N_gap_x_ + Nx_ + N_gap_x_));
  const double xmin_asmbly = r_low_.x() + static_cast<double>(index_assmbly_x) * assembly_dx_;
  if ((index_bin_x + 1) == N_gap_x_){
    return xmin_asmbly; // first portion of inter-assembly-gap inside the assembly
  } else if ((index_bin_x+1) <= (N_gap_x_ + Nx_)){ 
    // inside the fuel-pins
    return xmin_asmbly + 0.5 * inter_asmbly_gap_x_ * static_cast<double>(N_gap_x_) + static_cast<double>(index_bin_x-N_gap_x_) * dx_bin_; 
  } 

  // last portion of inter-assembly-gap inside the assembly
  return xmin_asmbly + assembly_dx_ - 0.5 * inter_asmbly_gap_x_ * static_cast<double>(N_gap_x_);
}

double RectAssemblyPositionFilter::x_max(const StaticVector3& index) const {
  if (real_Nx_ == 1) {
    return r_high_.x();
  }
  const std::size_t index_assmbly_x = static_cast<std::size_t>(std::floor(index[x_index_] / (N_gap_x_ + Nx_ + N_gap_x_)));
  const std::size_t index_bin_x = static_cast<std::size_t>(index[x_index_] - index_assmbly_x * (N_gap_x_ + Nx_ + N_gap_x_));
  const double xmin_asmbly = r_low_.x() + static_cast<double>(index_assmbly_x) * assembly_dx_;
  if ((index_bin_x + 1) == N_gap_x_){
    return xmin_asmbly + 0.5 * inter_asmbly_gap_x_ * static_cast<double>(N_gap_x_); // first portion of inter-assembly-gap inside the assembly
  } else if ((index_bin_x+1) <= (N_gap_x_ + Nx_)){ 
    // inside the fuel-pins
    return xmin_asmbly + 0.5 * inter_asmbly_gap_x_ * static_cast<double>(N_gap_x_) + static_cast<double>(index_bin_x-N_gap_x_) * dx_bin_ + dx_bin_; 
  } 

  // last portion of inter-assembly-gap inside the assembly
  return xmin_asmbly + assembly_dx_;
}

double RectAssemblyPositionFilter::y_min(const StaticVector3& index) const {
  if (real_Ny_ == 1) {
    return r_low_.y();
  }

  const std::size_t index_assmbly_y = static_cast<std::size_t>(std::floor(index[y_index_] / (N_gap_y_ + Ny_ + N_gap_y_)));
  const std::size_t index_bin_y = static_cast<std::size_t>(index[y_index_] - index_assmbly_y * (N_gap_y_ + Ny_ + N_gap_y_));
  const double ymin_asmbly = r_low_.y() + static_cast<double>(index_assmbly_y) * assembly_dy_;
  if ((index_bin_y + 1) == N_gap_y_){
    return ymin_asmbly; // first portion of inter-assembly-gap inside the assembly
  } else if ((index_bin_y+1) <= (N_gap_y_ + Ny_)){ 
    // inside the fuel-pins
    return ymin_asmbly + 0.5 * inter_asmbly_gap_y_ * static_cast<double>(N_gap_y_) + static_cast<double>(index_bin_y-N_gap_y_) * dy_bin_; 
  } 

  // last portion of inter-assembly-gap inside the assembly
  return ymin_asmbly + assembly_dy_ - 0.5 * inter_asmbly_gap_y_ * static_cast<double>(N_gap_y_);
}

double RectAssemblyPositionFilter::y_max(const StaticVector3& index) const {
  if (real_Ny_ == 1) return r_high_.y();
  
  const std::size_t index_assmbly_y = static_cast<std::size_t>(std::floor(index[y_index_] / (N_gap_y_ + Ny_ + N_gap_y_)));
  const std::size_t index_bin_y = static_cast<std::size_t>(index[y_index_] - index_assmbly_y * (N_gap_y_ + Ny_ + N_gap_y_));
  const double ymin_asmbly = r_low_.y() + static_cast<double>(index_assmbly_y) * assembly_dy_;
  if ((index_bin_y + 1) == N_gap_y_){
    return ymin_asmbly + 0.5 * inter_asmbly_gap_y_ * static_cast<double>(N_gap_y_); // first portion of inter-assembly-gap inside the assembly
  } else if ((index_bin_y+1) <= (N_gap_y_ + Ny_)){ 
    // inside the fuel-pins
    return ymin_asmbly + 0.5 * inter_asmbly_gap_y_ * static_cast<double>(N_gap_y_) + static_cast<double>(index_bin_y-N_gap_y_) * dy_bin_ + dy_bin_; 
  } 

  // last portion of inter-assembly-gap inside the assembly
  return ymin_asmbly + assembly_dy_;
}


double RectAssemblyPositionFilter::z_min(const StaticVector3& index) const {
  if (real_Nz_ == 1) return r_low_.z();
 
  const std::size_t index_assmbly_z = static_cast<std::size_t>(std::floor(index[z_index_] / (N_gap_z_ + Nz_ + N_gap_z_)));
  const std::size_t index_bin_z = static_cast<std::size_t>(index[z_index_] - index_assmbly_z * (N_gap_z_ + Nz_ + N_gap_z_));
  const double zmin_asmbly = r_low_.z() + static_cast<double>(index_assmbly_z) * assembly_dz_;
  if ((index_bin_z + 1) == N_gap_z_){
    return zmin_asmbly; // first portion of inter-assembly-gap inside the assembly
  } else if ((index_bin_z+1) <= (N_gap_z_ + Nz_)){ 
    // inside the fuel-pins
    return zmin_asmbly + 0.5 * inter_asmbly_gap_z_ * static_cast<double>(N_gap_z_) + static_cast<double>(index_bin_z-N_gap_z_) * dz_bin_; 
  } 

  // last portion of inter-assembly-gap inside the assembly
  return zmin_asmbly + assembly_dz_ - 0.5 * inter_asmbly_gap_z_ * static_cast<double>(N_gap_z_);
}

double RectAssemblyPositionFilter::z_max(const StaticVector3& index) const {
  if (real_Nz_ == 1) return r_high_.z();

  const std::size_t index_assmbly_z = static_cast<std::size_t>(std::floor(index[z_index_] / (N_gap_z_ + Nz_ + N_gap_z_)));
  const std::size_t index_bin_z = static_cast<std::size_t>(index[z_index_] - index_assmbly_z * (N_gap_z_ + Nz_ + N_gap_z_));
  const double zmin_asmbly = r_low_.z() + static_cast<double>(index_assmbly_z) * assembly_dz_;
  if ((index_bin_z + 1) == N_gap_z_){
    return zmin_asmbly + inter_asmbly_gap_z_ * static_cast<double>(N_gap_z_); // first portion of inter-assembly-gap inside the assembly
  } else if ((index_bin_z+1) <= (N_gap_z_ + Nz_)){ 
    // inside the fuel-pins
    return zmin_asmbly + 0.5 * inter_asmbly_gap_z_ * static_cast<double>(N_gap_z_) + static_cast<double>(index_bin_z-N_gap_z_) * dz_bin_ + dz_bin_; 
  } 

  // last portion of inter-assembly-gap inside the assembly
  return zmin_asmbly + assembly_dz_;
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
  grp.createAttribute("type", "rect-assembly-position-filter");

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
  std::array<std::size_t, 3> real_shape{real_Nx_, real_Ny_, real_Nz_};
  if (grp.hasAttribute("shape")) {
    grp.deleteAttribute("shape");
  }
  grp.createAttribute("shape", real_shape);

  // Save the shape of the assembly
  std::array<std::size_t, 3> assembly_shape{asmbly_Nx_, asmbly_Ny_, asmbly_Nz_};
  if (grp.hasAttribute("assembly-shape")){
    grp.deleteAttribute("assembly-shape");
  }
  grp.createAttribute("assembly-shape", assembly_shape);

  // save the inter assebly gap
  std::array<double, 3> inter_asmbly_gap{inter_asmbly_gap_x_, inter_asmbly_gap_y_, inter_asmbly_gap_z_};
  if (grp.hasAttribute("inter-assembly-gap")){
    grp.deleteAttribute("inter-assembly-gap");
  } 
  grp.createAttribute("inter-assembly-gap", inter_asmbly_gap);

  // save the shape of the bins excluding the gaps
  std::array<std::size_t, 3> bin_shape_in_asmbly{Nx_, Ny_, Nz_};
  if (grp.hasAttribute("bin-shape-per-assembly")){
    grp.deleteAttribute("bin-shape-per-assembly");
  }
  grp.createAttribute("bin-shape-per-assembly", bin_shape_in_asmbly);

  std::vector<double> x_bounds(real_Nx_ + 1, 0.);
  std::size_t itr = 0;
  double x0 = r_low_.x();
  x_bounds[itr] = x0; itr++;
  for (std::size_t i = 0; i < asmbly_Nx_; i++ ){
    if (N_gap_x_ == 1){
      x0 += 0.5 * inter_asmbly_gap_x_;
      x_bounds[itr] = x0; itr++;
    }

    for (std::size_t j = 0; j < Nx_; j++){
      x0 += dx_bin_;
      x_bounds[itr] = x0; itr++;
    }

    if (N_gap_x_ == 1){
      x0 += 0.5 * inter_asmbly_gap_x_;
      x_bounds[itr] = x0; itr++;
    }
  }
  grp.createDataSet("x-bounds", x_bounds);

  std::vector<double> y_bounds(real_Ny_ + 1, 0.);
  itr = 0;
  double y0 = r_low_.y();
  y_bounds[itr] = y0; itr++;
  for (std::size_t i = 0; i < asmbly_Ny_; i++ ){
    if (N_gap_y_ == 1){
      y0 += 0.5 * inter_asmbly_gap_y_;
      y_bounds[itr] = y0; itr++;
    }

    for (std::size_t j = 0; j < Ny_; j++){
      y0 += dy_bin_;
      y_bounds[itr] = y0; itr++;
    }

    if (N_gap_y_ == 1){
      y0 += 0.5 * inter_asmbly_gap_y_;
      y_bounds[itr] = y0; itr++;
    }
  }
  grp.createDataSet("y-bounds", y_bounds);

  std::vector<double> z_bounds(real_Nz_ + 1, 0.);
  itr = 0;
  double z0 = r_low_.z();
  z_bounds[itr] = z0; itr++;
  for (std::size_t i = 0; i < asmbly_Nz_; i++ ){
    if (N_gap_z_ == 1){
      z0 += 0.5 * inter_asmbly_gap_z_;
      z_bounds[itr] = z0; itr++;
    }

    for (std::size_t j = 0; j < Nz_; j++){
      z0 += dz_bin_;
      z_bounds[itr] = z0; itr++;
    }

    if (N_gap_z_ == 1){
      z0 += 0.5 * inter_asmbly_gap_z_;
      z_bounds[itr] = z0; itr++;
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
  } else if (!node["inter-assembly-gap"].IsSequence() || node["inter-assembly-gap"].size() != 3){
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
