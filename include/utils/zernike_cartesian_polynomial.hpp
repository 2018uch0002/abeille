#ifndef ZERNIKE_CARTESIAN_POLYNOMIAL
#define ZERNIKE_CARTESIAN_POLYNOMIAL

#include <utils/polynomial_expression.hpp>

//===================================
// Zernike Polynomials from order 0 up to order
// represented in cartesian coordinates
class ZernikeCartesianPolynomial {
 public:
  ZernikeCartesianPolynomial(const std::size_t order);

  // evaluation of all zernike-polynomial at distance t in direction (omega_x,
  // omega_y) from starting point(x0, y0)
  std::vector<double> evaluate_zernikes(const double t, const double omega_x,
                                        const double omega_y, const double x0,
                                        const double y0) const;

  // evaluation of all zernike polynomials at a point(x0, y0)
  std::vector<double> evaluate_zernikes(const double x0, const double y0) const;

  // evaluate the line-integration of all zernike polynomial
  std::vector<double> line_integrate_zernike(const double t_a, const double t_b,
                                             const double omega_x,
                                             const double omega_y,
                                             const double x0,
                                             const double y0) const;

 private:
  std::vector<PolynomialExpression> zr_poly_expression_;
  std::size_t order_, max_n_;

  // function for the factorial
  double factorial(std::size_t N) const {
    if (N == 1 || N == 0)
      return 1.;
    else
      return static_cast<double>(N) * factorial(N - 1);
  }

  // from the order get the n and l.
  // order = 0.5 * ( n*(n+2) + l )
  // assume the n = 0 and get the l, if not satisfied then increase the l
  // there will be a unique pair in theory
  std::pair<std::size_t, int> get_n_and_l(const std::size_t& order) const {
    int n = 0;
    while (n <= static_cast<int>(order)) {
      int l = 2 * static_cast<int>(order) - n * (n + 2);
      if (std::abs(l) <= n) return {static_cast<std::size_t>(n), l};
      n++;
    }
    return {0, 0};
  }
};

#endif