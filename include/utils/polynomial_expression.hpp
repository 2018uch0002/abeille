#ifndef POLYNOMIAL_EXPRESSION_H
#define POLYNOMIAL_EXPRESSION_H

#include <utils/error.hpp>

#include <iostream>
#include <string>
#include <vector>

class PolynomialExpression {
 public:
  enum class PolynomialSymbol { Dir_X, X0, Dir_Y, Y0, t };

  ~PolynomialExpression() = default;
  PolynomialExpression() = default;

  // monomial f(t) = coeff * t^pow;
  PolynomialExpression(double coeff, unsigned int pow_t);

  // linear polynomial along for x or y
  PolynomialExpression(std::string axis, double coeff = 1.);

  // construct the polynomial with the help of coefficient
  PolynomialExpression(const std::vector<double>& coeff,
                       const std::vector<unsigned int>& pow_t,
                       const std::vector<unsigned int>& pow_omega_x = {},
                       const std::vector<unsigned int>& pow_omega_y = {},
                       const std::vector<unsigned int>& pow_x0 = {},
                       const std::vector<unsigned int>& pow_y0 = {});

  PolynomialExpression(const PolynomialExpression&) = default;
  PolynomialExpression(PolynomialExpression&&) = default;

  // assigment operators
  PolynomialExpression& operator=(const PolynomialExpression&) = default;
  PolynomialExpression& operator=(PolynomialExpression&&) = default;

  // arithmetic operators
  PolynomialExpression operator+(const PolynomialExpression& poly_ex) const;
  void operator+=(const PolynomialExpression& poly_ex);
  PolynomialExpression operator-(const PolynomialExpression& poly_ex) const;
  void operator-=(const PolynomialExpression& poly_ex);
  PolynomialExpression operator*(const PolynomialExpression& poly_ex) const;
  void operator*=(const PolynomialExpression& poly_ex);
  PolynomialExpression operator*(const double val) const;
  void operator*=(const double val);

  // evalaute the function
  inline double evaluate(const double t, const double omega_x,
                         const double omega_y, const double x0,
                         const double y0) const {
    // first get the powers its max for t, omega_x, omega_y, x0, y0,
    std::vector<double> t_powers;
    t_powers.reserve(max_pow_t_ + 1);
    double t_m = 1.;
    for (std::size_t m = 0; m <= max_pow_t_; m++) {
      t_powers.push_back(t_m);
      t_m *= t;
    }

    std::vector<double> omega_x_powers;
    omega_x_powers.reserve(max_pow_omega_x_ + 1);
    double omega_x_m = 1.;
    for (std::size_t m = 0; m <= max_pow_omega_x_; m++) {
      omega_x_powers.push_back(omega_x_m);
      omega_x_m *= omega_x;
    }

    std::vector<double> omega_y_powers;
    omega_y_powers.reserve(max_pow_omega_y_ + 1);
    double omega_y_m = 1.;
    for (std::size_t m = 0; m <= max_pow_omega_y_; m++) {
      omega_y_powers.push_back(omega_y_m);
      omega_y_m *= omega_y;
    }

    std::vector<double> x0_powers;
    x0_powers.reserve(max_pow_x0_ + 1);
    double x0_m = 1.;
    for (std::size_t m = 0; m <= max_pow_x0_; m++) {
      x0_powers.push_back(x0_m);
      x0_m *= x0;
    }

    std::vector<double> y0_powers;
    y0_powers.reserve(max_pow_y0_ + 1);
    double y0_m = 1.;
    for (std::size_t m = 0; m <= max_pow_y0_; m++) {
      y0_powers.push_back(y0_m);
      y0_m *= y0;
    }

    double value = 0.;
    for (std::size_t i = 0; i < coeff_.size(); i++) {
      value += coeff_[i] * t_powers[pow_t_[i]] *
               omega_x_powers[pow_omega_x_[i]] *
               omega_y_powers[pow_omega_y_[i]] * x0_powers[pow_x0_[i]] *
               y0_powers[pow_y0_[i]];
    }
    return value;
  }

  // evaluate the function at (x0, y0), i.e., t = 0, omega_x = 0, omega_y = 0
  inline double evaluate_at_point(const double x0, const double y0) const {
    std::vector<double> x0_powers;
    x0_powers.reserve(max_pow_x0_ + 1);
    double x0_m = 1.;
    for (std::size_t m = 0; m <= max_pow_x0_; m++) {
      x0_powers.push_back(x0_m);
      x0_m *= x0;
    }

    std::vector<double> y0_powers;
    y0_powers.reserve(max_pow_y0_ + 1);
    double y0_m = 1.;
    for (std::size_t m = 0; m <= max_pow_y0_; m++) {
      y0_powers.push_back(y0_m);
      y0_m *= y0;
    }

    double value = 0.;
    for (std::size_t i = 0; i < max_loc_t_0_; i++) {
      value += coeff_[i] * x0_powers[pow_x0_[i]] * y0_powers[pow_y0_[i]];
    }
    return value;
  }

  // integrate the function w.r.t. to 't' from t_a to t_b
  inline double integrate(const double t_a, const double t_b,
                          const double omega_x, const double omega_y,
                          const double x0, const double y0) const {
    // first get the powers its max for t, omega_x, omega_y, x0, y0,
    std::vector<double> t_b_powers;
    t_b_powers.reserve(max_pow_t_ + 1);
    double t_b_m = t_b;
    for (std::size_t m = 0; m <= max_pow_t_; m++) {
      t_b_powers.push_back(t_b_m);
      t_b_m *= t_b;
    }

    std::vector<double> omega_x_powers;
    omega_x_powers.reserve(max_pow_omega_x_ + 1);
    double omega_x_m = 1.;
    for (std::size_t m = 0; m <= max_pow_omega_x_; m++) {
      omega_x_powers.push_back(omega_x_m);
      omega_x_m *= omega_x;
    }

    std::vector<double> omega_y_powers;
    omega_y_powers.reserve(max_pow_omega_y_ + 1);
    double omega_y_m = 1.;
    for (std::size_t m = 0; m <= max_pow_omega_y_; m++) {
      omega_y_powers.push_back(omega_y_m);
      omega_y_m *= omega_y;
    }

    std::vector<double> x0_powers;
    x0_powers.reserve(max_pow_x0_ + 1);
    double x0_m = 1.;
    for (std::size_t m = 0; m <= max_pow_x0_; m++) {
      x0_powers.push_back(x0_m);
      x0_m *= x0;
    }

    std::vector<double> y0_powers;
    y0_powers.reserve(max_pow_y0_ + 1);
    double y0_m = 1.;
    for (std::size_t m = 0; m <= max_pow_y0_; m++) {
      y0_powers.push_back(y0_m);
      y0_m *= y0;
    }

    double value = 0.;
    for (std::size_t i = 0; i < coeff_.size(); i++) {
      const double m_inv = 1. / static_cast<double>(pow_t_[i] + 1);
      value += m_inv * coeff_[i] * t_b_powers[pow_t_[i]] *
               omega_x_powers[pow_omega_x_[i]] *
               omega_y_powers[pow_omega_y_[i]] * x0_powers[pow_x0_[i]] *
               y0_powers[pow_y0_[i]];
    }

    if (t_a != 0.) {
      std::vector<double> t_a_powers;
      t_a_powers.reserve(max_pow_t_ + 1);
      double t_a_m = t_a;
      for (std::size_t m = 0; m <= max_pow_t_; m++) {
        t_a_powers.push_back(t_a_m);
        t_a_m *= t_a;
      }
      for (std::size_t i = 0; i < coeff_.size(); i++) {
        const double m_inv = 1. / static_cast<double>(pow_t_[i] + 1);
        value -= m_inv * coeff_[i] * t_a_powers[pow_t_[i]] *
                 omega_x_powers[pow_omega_x_[i]] *
                 omega_y_powers[pow_omega_y_[i]] * x0_powers[pow_x0_[i]] *
                 y0_powers[pow_y0_[i]];
      }
    }

    return value;
  }

  // to get the number of term
  std::size_t number_of_terms() const { return coeff_.size(); }

  // functions to get the different coeffs and pow
  const std::vector<double>& coefficients() const { return coeff_; }
  std::vector<double>& coefficients() { return coeff_; }
  const std::vector<unsigned int>& pow_t() const { return pow_t_; }
  std::vector<unsigned int>& pow_t() { return pow_t_; }
  const std::vector<unsigned int>& pow_omega_x() const { return pow_omega_x_; }
  std::vector<unsigned int>& pow_omega_x() { return pow_omega_x_; }
  const std::vector<unsigned int>& pow_omega_y() const { return pow_omega_y_; }
  std::vector<unsigned int>& pow_omega_y() { return pow_omega_y_; }
  const std::vector<unsigned int>& pow_x0() const { return pow_x0_; }
  std::vector<unsigned int>& pow_x0() { return pow_x0_; }
  const std::vector<unsigned int>& pow_y0() const { return pow_y0_; }
  std::vector<unsigned int>& pow_y0() { return pow_y0_; }

  // clear the polynomial
  void clear() {
    coeff_.clear();
    pow_t_.clear();
    pow_omega_x_.clear();
    pow_omega_y_.clear();
    pow_x0_.clear();
    pow_y0_.clear();
  }

  // print the expression
  void print_expression() const;

 private:
  std::vector<double> coeff_;
  std::vector<unsigned int> pow_t_;
  std::vector<unsigned int> pow_omega_x_;
  std::vector<unsigned int> pow_omega_y_;
  std::vector<unsigned int> pow_x0_;
  std::vector<unsigned int> pow_y0_;

  std::size_t max_pow_t_, max_pow_omega_x_, max_pow_omega_y_, max_pow_x0_,
      max_pow_y0_;
  std::size_t max_loc_t_0_;  // to determine upto how many terms have power to t
                             // is zero.

  // arranged into the assending order and avoid any repetion of terms
  public:
  void simplyfy_terms();
  void set_max_powers();  // pre-calculate max powers as it will help in
                          // function's evaluation
};

inline PolynomialExpression operator*(const double val,
                                      const PolynomialExpression& poly_ex) {
  return poly_ex * val;
}

#endif