#ifndef FUELPIN_APPROXIMATE_CANCELLATOR
#define FUELPIN_APPROXIMATE_CANCELLATOR

#include <cancelator/cancelator.hpp>

#include <unordered_map>

class FuelPinApproxCancelator : public Cancelator {
 public:
  enum class Orientation { X, Y, Z };  // to know the axial orientaion

  // If the number of rings are given.
  FuelPinApproxCancelator(Position low, Position high,
                          std::vector<double> assembly_pitch,
                          std::vector<std::uint32_t> assembly_shape,
                          std::vector<double> fuelpin_pitch,
                          std::vector<std::uint32_t> fuelpins_per_assembly,
                          std::vector<double> inter_assembly_gap, double radius,
                          std::uint32_t Nr, Orientation axial_orientation,
                          bool loop = false);

  // If the radius of the rings are given.
  FuelPinApproxCancelator(Position low, Position high,
                          std::vector<double> assembly_pitch,
                          std::vector<std::uint32_t> assembly_shape,
                          std::vector<double> fuelpin_pitch,
                          std::vector<std::uint32_t> fuelpins_per_assembly,
                          std::vector<double> inter_assembly_gap,
                          std::vector<double> radii,
                          Orientation axial_orientation, bool loop = false);

  // If the number of rings and the energy-groups are given.
  FuelPinApproxCancelator(Position low, Position high,
                          std::vector<double> assembly_pitch,
                          std::vector<std::uint32_t> assembly_shape,
                          std::vector<double> fuelpin_pitch,
                          std::vector<std::uint32_t> fuelpins_per_assembly,
                          std::vector<double> inter_assembly_gap, double radius,
                          std::uint32_t Nr, Orientation axial_orientation,
                          std::vector<double> energy_bounds, bool loop = false);

  // If the radius of the rings and the energy-groups are given.
  FuelPinApproxCancelator(Position low, Position high,
                          std::vector<double> assembly_pitch,
                          std::vector<std::uint32_t> assembly_shape,
                          std::vector<double> fuelpin_pitch,
                          std::vector<std::uint32_t> fuelpins_per_assembly,
                          std::vector<double> inter_assembly_gap,
                          std::vector<double> radii,
                          Orientation axial_orientation,
                          std::vector<double> energy_bounds, bool loop = false);

  bool add_particle(BankedParticle& p) override final;
  void perform_cancellation() override final;
  std::vector<BankedParticle> get_new_particles(RNG& /*rng*/) override final {
    return {};
  }
  void clear() override final { bins_.clear(); }
  void check_particle_mover_compatibility(
      const std::shared_ptr<IParticleMover>& /*pmover*/) const override final {}

  void write_output_info(H5::Group& grp) const override final;

 private:
  std::unordered_map<std::uint32_t, std::vector<BankedParticle*>> bins_;
  std::vector<double> energy_edges_;
  Position r_low_, r_high_;

  // assembly pitch
  double assembly_pitch_x_, assembly_pitch_y_, assembly_pitch_z_;
  double assembly_pitch_x_inv_, assembly_pitch_y_inv_, assembly_pitch_z_inv_;

  // fuel-pin pitch
  double pitch_x_, pitch_y_, pitch_z_;
  double pitch_x_inv_, pitch_y_inv_, pitch_z_inv_;

  // inter-assembly gaps
  double inter_asmbly_gap_x_, inter_asmbly_gap_y_, inter_asmbly_gap_z_;

  // squares of the radius of the rings
  std::vector<double> radius_ring_square_;

  // keep the number of fuel-pin in x, y, and z; then the number of rings, and
  // segments as well as energy-groups
  std::uint32_t Nx_asmbly_, Ny_asmbly_, Nz_asmbly_;
  std::uint32_t Nx_per_asmbly_, Ny_per_asmbly_, Nz_per_asmbly_;
  std::uint32_t Nr_, Ns_ = 4, Ne_ = 1;
  std::uint32_t Si_, Sj_, Sk_, Sr_, Ss_, Se_;  // Strides for indexing

  Orientation length_axis_ = Orientation::Z;

  bool loop;

  std::vector<std::uint32_t> sync_keys();
  void perform_cancellation_loop();
  // void perform_cancellation_vector();
  void perform_cancellation_full_vector();
};

std::shared_ptr<FuelPinApproxCancelator> make_fuelpin_approximate_cancelator(
    const YAML::Node& node);

#endif