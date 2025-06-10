#include <tallies/itally.hpp>
#include <tallies/lagrange_quad_element_fet.hpp>
#include <tallies/tallies.hpp>
#include <utils/output.hpp>

#include <boost/container/static_vector.hpp>
using StaticVector4 = boost::container::static_vector<size_t, 4>;

LagrangeQuadElementFET::LagrangeQuadElementFET(
    std::shared_ptr<CartesianFilter> position_filter,
    std::shared_ptr<EnergyFilter> energy_in, std::size_t polynomial_order,
    SpacialDomain sd, Quantity quantity, Estimator estimator, std::string name)
    : ITally(quantity, estimator, name),
      cartesian_filter_(position_filter),
      energy_in_(energy_in),
      poly_order_(polynomial_order),
      sd_(sd),
      index_x_(),
      index_y_(),
      index_z_(),
      loc_e_(0) {
  StaticVector4 tally_shape;
  // add the dimension for the energy_in_ only if exist
  if (energy_in_) {
    std::size_t ne = energy_in_->size();
    tally_shape.push_back(ne);
    loc_e_ = 1;
  }

  // get the shape or dimensions for the cartesian_filter_
  // the fet requires the co-ordinate information
  if (cartesian_filter_ == nullptr) {
    fatal_error("LagrangeQuadElementFET has nullptr cartesian-filter.");
  }

  // currently only 2D shape is supported.
  if (sd_ == SpacialDomain::XYZ) {
    fatal_error("LagrangeQuadElementFET only spports on 2D.");
  }

  StaticVector3 position_shape = cartesian_filter_->get_true_shape();

  index_x_ = 0;
  index_y_ = 1;
  index_z_ = 2;
  if (position_shape[0] == 1) {
    index_x_ = 0;
    index_y_--;
    index_z_--;
  }

  if (position_shape[1] == 1) {
    index_y_ = 0;
    index_z_--;
  }

  if (position_shape[2] == 1) {
    index_z_ = 0;
  }

  if (sd_ == SpacialDomain::XY) {
    tally_shape.push_back(position_shape[0] + 1);
    tally_shape.push_back(position_shape[1] + 1);
  } else {
    fatal_error("Only xy plane is supported.");
  }
  //   } else if (sd_ == SpacialDomain::YZ) {
  //     tally_shape.push_back(position_shape[1]);
  //     tally_shape.push_back(position_shape[2]);
  //   }
  //   if (sd_ == SpacialDomain::XZ) {
  //     tally_shape.push_back(position_shape[0]);
  //     tally_shape.push_back(position_shape[2]);
  //   }

  // currently only Lagrange Linear is supported
  if (poly_order_ != 1) {
    fatal_error("LagrangeQuadElementFET only spports linear shape function.");
  }

  // reallocate and fill with zeros for the tally avg, gen-score and variance
  tally_avg_.resize(tally_shape);
  tally_avg_.fill(0.0);

  tally_gen_score_.resize(tally_shape);
  tally_gen_score_.fill(0.0);

  tally_var_.resize(tally_shape);
  tally_var_.fill(0.0);
}

void LagrangeQuadElementFET::score_collision(const Particle& p,
                                             const Tracker& trkr,
                                             MaterialHelper& mat) {
  StaticVector4 indices;
  // get the energy-index, if energy-filter exists
  if (energy_in_) {
    std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());
    if (E_indx.has_value() == false) {
      // Not inside any energy bin. Don't score.
      return;
    }

    indices.push_back(E_indx.value());
  }

  // get the cartisian_filter indices
  StaticVector3 position_index = cartesian_filter_->get_indices(trkr);
  if (position_index.empty()) {
    // No bin is found, don't score.
    return;
  }

  indices.push_back(position_index[index_x_]);
  indices.push_back(position_index[index_y_]);

  const double Et = mat.Et(p.E());
  const double collision_score =
      particle_base_score(p.E(), p.wgt(), p.wgt2(), &mat) / Et;

  // add the score at the xmin-ymin
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  tally_gen_score_.element(indices.begin(), indices.end()) += collision_score;

  // add the score at the xmax-ymin
  indices[loc_e_] += 1;
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  tally_gen_score_.element(indices.begin(), indices.end()) += collision_score;

  // add the score at the xmax-ymax
  indices[loc_e_ + 1] += 1;
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  tally_gen_score_.element(indices.begin(), indices.end()) += collision_score;

  // add the score at the xmin-ymax
  indices[loc_e_] -= 1;
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  tally_gen_score_.element(indices.begin(), indices.end()) += collision_score;

std::cout << "Check wheather we correctly tallying or not with collision socore = " << collision_score << std::endl;
for (auto&p : position_index)
  std::cout << p << "\t";
std::cout << "\n"<< std::endl;

std::vector<std::size_t> itr(tally_gen_score_.shape().begin(), tally_gen_score_.shape().end());
for (auto&p : itr)
  std::cout << p << "\t";
std::cout << "\n"<< std::endl;

for(std::size_t ix : {0, 1, 2, 3}){
  for(std::size_t iy : {0, 1, 2, 3}){
    std::cout << "ix = " << ix << ", iy = " << iy << ":\t" << tally_gen_score_(ix, iy) << std::endl;
  }
  std::cout << "---------" << std::endl;
}
fatal_error("JOB DONE!");
}

std::string LagrangeQuadElementFET::spacial_domain() const {
  switch (sd_) {
    case SpacialDomain::XY:
      return "xy";
      break;

    case SpacialDomain::YZ:
      return "yz";
      break;

    case SpacialDomain::XZ:
      return "xz";
      break;

    case SpacialDomain::XYZ:
      return "xyz";
      break;

    default:
      return "unknown";
  }
}

void LagrangeQuadElementFET::write_tally() {
  // Only master can write tallies, as only master has a copy
  // of the mean and variance.
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  // Create the group for the tally
  auto tally_grp = h5.createGroup("results/" + tally_name_);

  // Save the type
  tally_grp.createAttribute("type", "lagrange-quad-element-fet");

  // Save the quantity
  tally_grp.createAttribute("quantity", quantity_str());
  if (quantity_.type == Quantity::Type::MT) {
    tally_grp.createAttribute("mt", quantity_.mt);
  }

  // Save the polynomial-order
  tally_grp.createAttribute("polynomial-order", poly_order_);

  // Save the spatial-domain
  tally_grp.createAttribute("spatial-domain", spacial_domain());

  // Save the estimator
  tally_grp.createAttribute("estimator", estimator_str());

  // Save energy-in filter id
  if (energy_in_) {
    tally_grp.createAttribute("energy-filter", energy_in_->id());
  }

  // Save position filter id
  tally_grp.createAttribute("position-filter", cartesian_filter_->id());

  // Convert flux_var to the error on the mean
  this->var_to_std_on_mean();

  // Add data sets for the average and the standard deviation
  std::vector<std::size_t> shape(tally_avg_.shape().begin(),
                                 tally_avg_.shape().end());
  auto avg_dset = tally_grp.createDataSet<double>("avg", H5::DataSpace(shape));
  avg_dset.write_raw(tally_avg_.data());

  auto std_dset = tally_grp.createDataSet<double>("std", H5::DataSpace(shape));
  std_dset.write_raw(tally_var_.data());
}

std::shared_ptr<LagrangeQuadElementFET> make_lagrange_quad_element_fet(
    const YAML::Node& node) {
  // Check the name of the tally is given or not.
  if (!node["name"] || node["name"].IsScalar() == false) {
    fatal_error("No valid name is provided on tally.");
  }
  std::string name = node["name"].as<std::string>();

  // Check for the quantity
  std::string given_quantity = "";
  if (!node["quantity"] || node["quantity"].IsScalar() == false) {
    fatal_error("No valid quantity entry is given for tally " + name + ".");
  }
  given_quantity = node["quantity"].as<std::string>();
  Quantity quant = read_quantity(node, name);
  const bool source_like = (quant.type == Quantity::Type::Source ||
                            quant.type == Quantity::Type::RealSource ||
                            quant.type == Quantity::Type::ImagSource);

  std::string estimator_name = "collision";
  if (source_like || estimator_name != "collision") {
    fatal_error("No valid estimator entry is given for tally " + name + ".");
  }

  // Check the estimator is given or not. We default to collision estimators.
  if (node["estimator"] && node["estimator"].IsScalar()) {
    estimator_name = node["estimator"].as<std::string>();
  } else if (node["estimator"]) {
    fatal_error("Invalid estimator entry is given on tally " + name + ".");
  }

  Estimator estimator;
  if (estimator_name == "collision") {
    estimator = Estimator::Collision;
  } else if (estimator_name == "track-length") {
    estimator = Estimator::TrackLength;
  } else if (estimator_name == "source") {
    estimator = Estimator::Source;
  } else {
    fatal_error("Invalid estimator is given for tally " + name + ".");
  }

  if (source_like && estimator != Estimator::Source) {
    std::stringstream mssg;
    mssg << "Tally " << name
         << " has a source-like quantity but not a source estimator.";
    fatal_error(mssg.str());
  }

  // get the polynomaial order
  std::size_t poly_order = 1;
  if (node["polynomial-order"] && node["polynomial-order"].IsScalar()) {
    poly_order = node["polynomial-order"].as<std::size_t>();
    if (poly_order != 1) {
      std::stringstream mssg;
      mssg << "Tally " << name
           << " current only supports the linear lagrange shape element.";
      fatal_error(mssg.str());
    }
  } else if (node["polynomial-order"]) {
    std::stringstream mssg;
    mssg << "Tally " << name << " has invalid polynomial order.";
    fatal_error(mssg.str());
  }

  // get the spatial domain
  LagrangeQuadElementFET::SpacialDomain sd =
      LagrangeQuadElementFET::SpacialDomain::XY;
  if (node["spatial-domain"] && node["spatial-domain"].IsScalar()) {
    std::string sd_name = node["spatial-domain"].as<std::string>();
    if (sd_name != "xy") {
      std::stringstream mssg;
      mssg << "Tally " << name << " only supports xy spatial-domain.";
      fatal_error(mssg.str());
    }
  } else if (node["spatial-domain"]) {
    std::stringstream mssg;
    mssg << "Tally " << name << " has invalid spatial-domain.";
    fatal_error(mssg.str());
  }

  // get the tallies instance
  auto& tallies = Tallies::instance();

  // Get the energy bounds, if any is given
  std::shared_ptr<EnergyFilter> energy_filter = nullptr;
  if (node["energy-filter"] && node["energy-filter"].IsScalar()) {
    std::size_t energy_id = node["energy-filter"].as<std::size_t>();
    energy_filter = tallies.get_energy_filter(energy_id);
    if (energy_filter == nullptr) {
      std::stringstream mssg;
      mssg << "For tally " << name << ", cannot find energy filter with id "
           << energy_id << ".";
      fatal_error(mssg.str());
    }
  } else if (node["energy-filter"]) {
    fatal_error("Invalid energy-filter entry on tally " + name + ".");
  }

  // Get the position filter
  std::shared_ptr<CartesianFilter> position_filter = nullptr;
  if (node["position-filter"] && node["position-filter"].IsScalar()) {
    std::size_t position_id = node["position-filter"].as<std::size_t>();
    position_filter = tallies.get_cartesian_filter(position_id);
    if (position_filter == nullptr) {
      std::stringstream mssg;
      mssg << "For tally " << name << ", cannot find position filter with id "
           << position_id << ".";
      fatal_error(mssg.str());
    }
  }

  // For the general tally
  std::shared_ptr<LagrangeQuadElementFET> tally =
      std::make_shared<LagrangeQuadElementFET>(position_filter, energy_filter,
                                               poly_order, sd, quant, estimator,
                                               name);

  return tally;
}