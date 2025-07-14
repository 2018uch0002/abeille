#ifndef GAUSS_LEGENDRE_QUADRATURE_H
#define GAUSS_LEGENDRE_QUADRATURE_H

#include <array>
#include <span>
#include <string>

struct Abscissa_and_Weight {
  double abscissa = 0.;
  double weight = 0.;
};

template <std::size_t N>
class GaussLegendreQuad {
 public:
  std::span<const Abscissa_and_Weight> quadrature_set() const {
    return {abcisa_wgt_.begin(), abcisa_wgt_.end()};
  }

  std::string type_str() const { return "gauss-legendre"; }

 private:
  static const std::array<Abscissa_and_Weight, N> abcisa_wgt_;
};

template <std::size_t N>
class MidPointQuad {
 public:
  std::span<const Abscissa_and_Weight> quadrature_set() const {
    return {abcisa_wgt_.begin(), abcisa_wgt_.end()};
  }

  std::string type_str() const { return "mid-point"; }

 private:
  static const std::array<Abscissa_and_Weight, N> abcisa_wgt_;
};

#endif
