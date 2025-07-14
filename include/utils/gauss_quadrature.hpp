#ifndef GAUSS_QUADRATURE
#define GAUSS_QUADRATURE

#include <utils/gauss_legendre_quadrature.hpp>

#include <variant>

using QuadratureType = std::variant<
    GaussLegendreQuad<1>, GaussLegendreQuad<2>, GaussLegendreQuad<3>,
    GaussLegendreQuad<4>, GaussLegendreQuad<5>, GaussLegendreQuad<6>,
    GaussLegendreQuad<7>, GaussLegendreQuad<8>, GaussLegendreQuad<9>,
    GaussLegendreQuad<10>, GaussLegendreQuad<11>, GaussLegendreQuad<12>,
    GaussLegendreQuad<16>, GaussLegendreQuad<20>, GaussLegendreQuad<32>,
    GaussLegendreQuad<64>, MidPointQuad<10>>;

class GaussQuadrature {
 public:
  GaussQuadrature(QuadratureType gq) : gq_(gq), abcisa_wgt_() {
    abcisa_wgt_ =
        std::visit([](const auto& gq) { return gq.quadrature_set(); }, gq_);
  }

  const std::span<const Abscissa_and_Weight>& quadrature_set() const {
    return abcisa_wgt_;
  }

  std::size_t size() const { return abcisa_wgt_.size(); }
  std::string type_str() const {
    std::string type_str_ =
        std::visit([](const auto& gq) { return gq.type_str(); }, gq_);
    return type_str_;
  }

 private:
  QuadratureType gq_;
  std::span<const Abscissa_and_Weight> abcisa_wgt_;
};

#endif