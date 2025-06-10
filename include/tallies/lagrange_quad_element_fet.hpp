#ifndef LAGRANGE_QUAD_ELEMENT_FET_H
#define LAGRANGE_QUAD_ELEMENT_FET_H

#include <tallies/cartesian_filter.hpp>
#include <tallies/energy_filter.hpp>
#include <tallies/itally.hpp>

#include <boost/container/static_vector.hpp>
using StaticVector4 = boost::container::static_vector<size_t, 4>;

class LagrangeQuadElementFET : public ITally {
 public:
  enum class SpacialDomain { XY, YZ, XZ, XYZ };

  LagrangeQuadElementFET(std::shared_ptr<CartesianFilter> position_filter,
                         std::shared_ptr<EnergyFilter> energy_in,
                         std::size_t polynomial_order, SpacialDomain sd,
                         Quantity quantity, Estimator estimator,
                         std::string name);

  ~LagrangeQuadElementFET() = default;

  void score_collision(const Particle& p, const Tracker& trkr,
                       MaterialHelper& mat) override final;

  void score_flight(const Particle& /*p*/, const Tracker& /*trkr*/,
                    double /*d_flight*/,
                    MaterialHelper& /*mat*/) override final {
    fatal_error("the track-length for the legendre-fet is not supoorted yet.");
  }

  void score_source(const BankedParticle& /*p*/) override final {
    fatal_error("the track-length for the legendre-fet is not supoorted yet.");
  }

  double evaluate(const Position& /*r*/, const double& /*E*/) const override final {
    fatal_error("the track-length for the legendre-fet is not supoorted yet.");
    return {};
  }
  std::vector<double> evaluate(
      const std::vector<std::pair<Position, double>> /*r_E*/) const override final {
    fatal_error("the track-length for the legendre-fet is not supoorted yet.");
    return {};
  }

  std::size_t get_polynomial_order() { return poly_order_; }

  std::string spacial_domain() const;

  void write_tally() override final;

 private:
  std::shared_ptr<CartesianFilter> cartesian_filter_;
  std::shared_ptr<EnergyFilter> energy_in_;
  std::size_t poly_order_ = 1;
  SpacialDomain sd_;

  std::size_t index_x_, index_y_, index_z_, loc_e_;
};

std::shared_ptr<LagrangeQuadElementFET> make_lagrange_quad_element_fet(
    const YAML::Node& node);

#endif