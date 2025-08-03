#include <utils/polynomial_expression.hpp>

#include <sstream>

PolynomialExpression::PolynomialExpression(double coeff, unsigned int pow_t)
    : coeff_({coeff}),
      pow_t_({pow_t}),
      pow_omega_x_({0}),
      pow_omega_y_({0}),
      pow_x0_({0}),
      pow_y0_({0}),
      max_pow_t_(),
      max_pow_omega_x_(),
      max_pow_omega_y_(),
      max_pow_x0_(),
      max_pow_y0_(),
      max_loc_t_0_() {
  simplyfy_terms();
  set_max_powers();
}

PolynomialExpression::PolynomialExpression(std::string axis, double coeff)
    : coeff_({coeff, coeff}),
      pow_t_({1, 0}),
      pow_omega_x_({0, 0}),
      pow_omega_y_({0, 0}),
      pow_x0_({0, 0}),
      pow_y0_({0, 0}),
      max_pow_t_(),
      max_pow_omega_x_(),
      max_pow_omega_y_(),
      max_pow_x0_(),
      max_pow_y0_(),
      max_loc_t_0_() {
  if (axis == "X") {
    pow_omega_x_[0] = 1;
    pow_x0_[1] = 1;
  } else if (axis == "Y") {
    pow_omega_y_[0] = 1;
    pow_y0_[1] = 1;
  } else {
    fatal_error("Invalid axis is given in PolynomialExpression.");
  }
  simplyfy_terms();
  set_max_powers();
}

PolynomialExpression::PolynomialExpression(
    const std::vector<double>& coeff, const std::vector<unsigned int>& pow_t,
    const std::vector<unsigned int>& pow_omega_x,
    const std::vector<unsigned int>& pow_omega_y,
    const std::vector<unsigned int>& pow_x0,
    const std::vector<unsigned int>& pow_y0)
    : coeff_(coeff),
      pow_t_(pow_t),
      pow_omega_x_(pow_omega_x),
      pow_omega_y_(pow_omega_y),
      pow_x0_(pow_x0),
      pow_y0_(pow_y0),
      max_pow_t_(),
      max_pow_omega_x_(),
      max_pow_omega_y_(),
      max_pow_x0_(),
      max_pow_y0_(),
      max_loc_t_0_() {
  simplyfy_terms();
  set_max_powers();
}

void PolynomialExpression::print_expression() const {
  std::stringstream expression;
  for (std::size_t i = 0; i < pow_t_.size(); i++) {
    if (coeff_[i] != 0.) {
      if (i != 0) {
        if (coeff_[i] >= 0.)
          expression << " + ";
        else
          expression << " ";
      }

      if (coeff_[i] != 1.) {
        expression << coeff_[i] << " * ";
      }

      if (pow_omega_x_[i] != 0)
        expression << "(omega_x ^ " << pow_omega_x_[i] << ")" << " * ";

      if (pow_omega_y_[i] != 0)
        expression << "(omega_y ^ " << pow_omega_y_[i] << ")" << " * ";

      if (pow_x0_[i] != 0) expression << "(x0 ^ " << pow_x0_[i] << ")" << " * ";

      if (pow_y0_[i] != 0) expression << "(y0 ^ " << pow_y0_[i] << ")" << " * ";

      expression << "(t ^ " << pow_t_[i] << ")";
    }
    expression << "\n";
  }
  std::cout << "Expression:\n " << expression.str() << std::endl;
}

PolynomialExpression PolynomialExpression::operator+(
    const PolynomialExpression& poly_ex) const {
  const std::size_t num_terms = coeff_.size() + poly_ex.number_of_terms();
  std::vector<double> new_coeff = coeff_;
  new_coeff.reserve(num_terms);
  new_coeff.insert(new_coeff.end(), poly_ex.coefficients().begin(),
                   poly_ex.coefficients().end());

  std::vector<unsigned int> new_pow_t = pow_t_;
  new_pow_t.reserve(num_terms);
  new_pow_t.insert(new_pow_t.end(), poly_ex.pow_t().begin(),
                   poly_ex.pow_t().end());

  std::vector<unsigned int> new_pow_omega_x = pow_omega_x_;
  new_pow_omega_x.reserve(num_terms);
  new_pow_omega_x.insert(new_pow_omega_x.end(), poly_ex.pow_omega_x().begin(),
                         poly_ex.pow_omega_x().end());

  std::vector<unsigned int> new_pow_omega_y = pow_omega_y_;
  new_pow_omega_y.reserve(num_terms);
  new_pow_omega_y.insert(new_pow_omega_y.end(), poly_ex.pow_omega_y().begin(),
                         poly_ex.pow_omega_y().end());

  std::vector<unsigned int> new_pow_x0 = pow_x0_;
  new_pow_x0.reserve(num_terms);
  new_pow_x0.insert(new_pow_x0.end(), poly_ex.pow_x0().begin(),
                    poly_ex.pow_x0().end());

  std::vector<unsigned int> new_pow_y0 = pow_y0_;
  new_pow_y0.reserve(num_terms);
  new_pow_y0.insert(new_pow_y0.end(), poly_ex.pow_y0().begin(),
                    poly_ex.pow_y0().end());

  return PolynomialExpression(new_coeff, new_pow_t, new_pow_omega_x,
                              new_pow_omega_y, new_pow_x0, new_pow_y0);
}

void PolynomialExpression::operator+=(const PolynomialExpression& poly_ex) {
  const std::size_t num_terms = coeff_.size() + poly_ex.number_of_terms();
  coeff_.reserve(num_terms);
  coeff_.insert(coeff_.end(), poly_ex.coefficients().begin(),
                poly_ex.coefficients().end());

  pow_t_.reserve(num_terms);
  pow_t_.insert(pow_t_.end(), poly_ex.pow_t().begin(), poly_ex.pow_t().end());

  pow_omega_x_.reserve(num_terms);
  pow_omega_x_.insert(pow_omega_x_.end(), poly_ex.pow_omega_x().begin(),
                      poly_ex.pow_omega_x().end());

  pow_omega_y_.reserve(num_terms);
  pow_omega_y_.insert(pow_omega_y_.end(), poly_ex.pow_omega_y().begin(),
                      poly_ex.pow_omega_y().end());

  pow_x0_.reserve(num_terms);
  pow_x0_.insert(pow_x0_.end(), poly_ex.pow_x0().begin(),
                 poly_ex.pow_x0().end());

  pow_y0_.reserve(num_terms);
  pow_y0_.insert(pow_y0_.end(), poly_ex.pow_y0().begin(),
                 poly_ex.pow_y0().end());

  simplyfy_terms();
  set_max_powers();
}

PolynomialExpression PolynomialExpression::operator-(
    const PolynomialExpression& poly_ex) const {
  const std::size_t num_terms = coeff_.size() + poly_ex.number_of_terms();
  std::vector<double> new_coeff = coeff_;
  new_coeff.reserve(num_terms);
  const std::size_t i0 = coeff_.size();
  new_coeff.insert(new_coeff.end(), poly_ex.coefficients().begin(),
                   poly_ex.coefficients().end());
  for (std::size_t i = i0; i < new_coeff.size(); i++) {
    new_coeff[i] *= -1.;
  }

  std::vector<unsigned int> new_pow_t = pow_t_;
  new_pow_t.reserve(num_terms);
  new_pow_t.insert(new_pow_t.end(), poly_ex.pow_t().begin(),
                   poly_ex.pow_t().end());

  std::vector<unsigned int> new_pow_omega_x = pow_omega_x_;
  new_pow_omega_x.reserve(num_terms);
  new_pow_omega_x.insert(new_pow_omega_x.end(), poly_ex.pow_omega_x().begin(),
                         poly_ex.pow_omega_x().end());

  std::vector<unsigned int> new_pow_omega_y = pow_omega_y_;
  new_pow_omega_y.reserve(num_terms);
  new_pow_omega_y.insert(new_pow_omega_y.end(), poly_ex.pow_omega_y().begin(),
                         poly_ex.pow_omega_y().end());

  std::vector<unsigned int> new_pow_x0 = pow_x0_;
  new_pow_x0.reserve(num_terms);
  new_pow_x0.insert(new_pow_x0.end(), poly_ex.pow_x0().begin(),
                    poly_ex.pow_x0().end());

  std::vector<unsigned int> new_pow_y0 = pow_y0_;
  new_pow_y0.reserve(num_terms);
  new_pow_y0.insert(new_pow_y0.end(), poly_ex.pow_y0().begin(),
                    poly_ex.pow_y0().end());

  return PolynomialExpression(new_coeff, new_pow_t, new_pow_omega_x,
                              new_pow_omega_y, new_pow_x0, new_pow_y0);
}

void PolynomialExpression::operator-=(const PolynomialExpression& poly_ex) {
  const std::size_t num_terms = coeff_.size() + poly_ex.number_of_terms();
  coeff_.reserve(num_terms);
  const std::size_t i0 = coeff_.size();
  coeff_.insert(coeff_.end(), poly_ex.coefficients().begin(),
                poly_ex.coefficients().end());
  for (std::size_t i = i0; i < coeff_.size(); i++) {
    coeff_[i] *= -1.;
  }

  pow_t_.reserve(num_terms);
  pow_t_.insert(pow_t_.end(), poly_ex.pow_t().begin(), poly_ex.pow_t().end());

  pow_omega_x_.reserve(num_terms);
  pow_omega_x_.insert(pow_omega_x_.end(), poly_ex.pow_omega_x().begin(),
                      poly_ex.pow_omega_x().end());

  pow_omega_y_.reserve(num_terms);
  pow_omega_y_.insert(pow_omega_y_.end(), poly_ex.pow_omega_y().begin(),
                      poly_ex.pow_omega_y().end());

  pow_x0_.reserve(num_terms);
  pow_x0_.insert(pow_x0_.end(), poly_ex.pow_x0().begin(),
                 poly_ex.pow_x0().end());

  pow_y0_.reserve(num_terms);
  pow_y0_.insert(pow_y0_.end(), poly_ex.pow_y0().begin(),
                 poly_ex.pow_y0().end());

  simplyfy_terms();
  set_max_powers();
}

PolynomialExpression PolynomialExpression::operator*(
    const PolynomialExpression& poly_ex) const {
  const std::size_t current_term = coeff_.size();
  const std::size_t poly_ex_terms = poly_ex.number_of_terms();
  const std::size_t num_terms = current_term * poly_ex.number_of_terms();

  // for new expression
  std::vector<double> new_coeff;
  new_coeff.reserve(num_terms);
  std::vector<unsigned int> new_pow_t;
  new_pow_t.reserve(num_terms);
  std::vector<unsigned int> new_pow_omega_x;
  new_pow_omega_x.reserve(num_terms);
  std::vector<unsigned int> new_pow_omega_y;
  new_pow_omega_y.reserve(num_terms);
  std::vector<unsigned int> new_pow_x0;
  new_pow_x0.reserve(num_terms);
  std::vector<unsigned int> new_pow_y0;
  new_pow_y0.reserve(num_terms);

  // first we will add the extra terms, that means starting from the second
  // term of poly_ex wil be push_back. When all the terms will be added,
  // then, the first term of poly_ex will be multiplied in place.
  for (std::size_t i = 0; i < poly_ex_terms; i++) {
    const double multiply_coeff = poly_ex.coefficients()[i];
    const unsigned int multiply_pow_t = poly_ex.pow_t()[i];
    const unsigned int multiply_pow_omega_x = poly_ex.pow_omega_x()[i];
    const unsigned int multiply_pow_omega_y = poly_ex.pow_omega_y()[i];
    const unsigned int multiply_pow_x0 = poly_ex.pow_x0()[i];
    const unsigned int multiply_pow_y0 = poly_ex.pow_y0()[i];

    for (std::size_t j = 0; j < current_term; j++) {
      new_coeff.push_back(coeff_[j] * multiply_coeff);
      new_pow_t.push_back(pow_t_[j] + multiply_pow_t);
      new_pow_omega_x.push_back(pow_omega_x_[j] + multiply_pow_omega_x);
      new_pow_omega_y.push_back(pow_omega_y_[j] + multiply_pow_omega_y);
      new_pow_x0.push_back(pow_x0_[j] + multiply_pow_x0);
      new_pow_y0.push_back(pow_y0_[j] + multiply_pow_y0);
    }
  }

  return PolynomialExpression(new_coeff, new_pow_t, new_pow_omega_x,
                              new_pow_omega_y, new_pow_x0, new_pow_y0);
}

void PolynomialExpression::operator*=(const PolynomialExpression& poly_ex) {
  const std::size_t current_term = coeff_.size();
  const std::size_t num_terms = current_term * poly_ex.number_of_terms();

  // reserve the memory
  coeff_.reserve(num_terms);
  pow_t_.reserve(num_terms);
  pow_omega_x_.reserve(num_terms);
  pow_omega_y_.reserve(num_terms);
  pow_x0_.reserve(num_terms);
  pow_y0_.reserve(num_terms);

  // first we will add the extra terms, that means starting from the second
  // term of poly_ex wil be push_back. When all the terms will be added,
  // then, the first term of poly_ex will be multiplied in place.
  for (std::size_t i = 1; i < poly_ex.number_of_terms(); i++) {
    const double multiply_coeff = poly_ex.coefficients()[i];
    const unsigned int multiply_pow_t = poly_ex.pow_t()[i];
    const unsigned int multiply_pow_omega_x = poly_ex.pow_omega_x()[i];
    const unsigned int multiply_pow_omega_y = poly_ex.pow_omega_y()[i];
    const unsigned int multiply_pow_x0 = poly_ex.pow_x0()[i];
    const unsigned int multiply_pow_y0 = poly_ex.pow_y0()[i];

    for (std::size_t j = 0; j < current_term; j++) {
      coeff_.push_back(coeff_[j] * multiply_coeff);
      pow_t_.push_back(pow_t_[j] + multiply_pow_t);
      pow_omega_x_.push_back(pow_omega_x_[j] + multiply_pow_omega_x);
      pow_omega_y_.push_back(pow_omega_y_[j] + multiply_pow_omega_y);
      pow_x0_.push_back(pow_x0_[j] + multiply_pow_x0);
      pow_y0_.push_back(pow_y0_[j] + multiply_pow_y0);
    }
  }

  // modifiy the coeff. and powers in place
  const double multiply_coeff = poly_ex.coefficients()[0];
  const unsigned int multiply_pow_t = poly_ex.pow_t()[0];
  const unsigned int multiply_pow_omega_x = poly_ex.pow_omega_x()[0];
  const unsigned int multiply_pow_omega_y = poly_ex.pow_omega_y()[0];
  const unsigned int multiply_pow_x0 = poly_ex.pow_x0()[0];
  const unsigned int multiply_pow_y0 = poly_ex.pow_y0()[0];
  for (std::size_t j = 0; j < current_term; j++) {
    coeff_[j] *= multiply_coeff;
    pow_t_[j] += multiply_pow_t;
    pow_omega_x_[j] += multiply_pow_omega_x;
    pow_omega_y_[j] += multiply_pow_omega_y;
    pow_x0_[j] += multiply_pow_x0;
    pow_y0_[j] += multiply_pow_y0;
  }

  simplyfy_terms();
  set_max_powers();
}

PolynomialExpression PolynomialExpression::operator*(const double val) const {
  std::vector<double> modify_coeff;
  modify_coeff.reserve(coeff_.size());
  for (std::size_t i = 0; i < coeff_.size(); i++) {
    modify_coeff.push_back(coeff_[i] * val);
  }
  return PolynomialExpression(modify_coeff, pow_t_, pow_omega_x_, pow_omega_y_,
                              pow_x0_, pow_y0_);
}

void PolynomialExpression::operator*=(const double val) {
  for (std::size_t i = 0; i < coeff_.size(); i++) {
    coeff_[i] *= val;
  }
}

void PolynomialExpression::simplyfy_terms() {
  // first sort the powers in ascending order for t and adjust the coefficents
  // if two terms are excalty same.

  std::size_t N_terms = coeff_.size();

  for (std::size_t i = 1; i < N_terms; i++) {
    const unsigned int key_pow_t = pow_t_[i];
    const unsigned int key_pow_omega_x = pow_omega_x_[i];
    const unsigned int key_pow_omega_y = pow_omega_y_[i];
    const unsigned int key_pow_x0 = pow_x0_[i];
    const unsigned int key_pow_y0 = pow_y0_[i];
    const double key_coeff = coeff_[i];

    std::size_t uj = i;
    while (uj > 0 && (key_pow_t < pow_t_[uj - 1])) {
      pow_t_[uj] = pow_t_[uj - 1];
      coeff_[uj] = coeff_[uj - 1];
      pow_omega_x_[uj] = pow_omega_x_[uj - 1];
      pow_omega_y_[uj] = pow_omega_y_[uj - 1];
      pow_x0_[uj] = pow_x0_[uj - 1];
      pow_y0_[uj] = pow_y0_[uj - 1];
      uj--;
    }
    pow_t_[uj] = key_pow_t;
    pow_omega_x_[uj] = key_pow_omega_x;
    pow_omega_y_[uj] = key_pow_omega_y;
    pow_x0_[uj] = key_pow_x0;
    pow_y0_[uj] = key_pow_y0;
    coeff_[uj] = key_coeff;

    if (uj > 0) {
      if (pow_t_[uj] != pow_t_[uj - 1]) continue;
      if (pow_omega_x_[uj] != pow_omega_x_[uj - 1]) continue;
      if (pow_omega_y_[uj] != pow_omega_y_[uj - 1]) continue;
      if (pow_x0_[uj] != pow_x0_[uj - 1]) continue;
      if (pow_y0_[uj] != pow_y0_[uj - 1]) continue;

      pow_t_.erase(pow_t_.begin() + static_cast<std::size_t>(uj));
      pow_omega_x_.erase(pow_omega_x_.begin() + static_cast<std::size_t>(uj));
      pow_omega_y_.erase(pow_omega_y_.begin() + static_cast<std::size_t>(uj));
      pow_x0_.erase(pow_x0_.begin() + static_cast<std::size_t>(uj));
      pow_y0_.erase(pow_y0_.begin() + static_cast<std::size_t>(uj));
      coeff_[uj - 1] += coeff_[uj];
      coeff_.erase(coeff_.begin() + static_cast<std::size_t>(uj));
      N_terms--;
      i--;
    }
  }

  // set for the location, where upto t^0
  max_loc_t_0_ = 0;
  for (std::size_t i = 0; i < pow_t_.size(); i++) {
    if (pow_t_[i] == 0) {
      max_loc_t_0_++;
    } else {
      break;
    }
  }
}

void PolynomialExpression::set_max_powers() {
  max_pow_t_ = pow_t_[0];
  for (std::size_t i = 1; i < pow_t_.size(); i++) {
    if (pow_t_[i] > max_pow_t_) max_pow_t_ = pow_t_[i];
  }

  max_pow_omega_x_ = pow_omega_x_[0];
  for (std::size_t i = 1; i < pow_omega_x_.size(); i++) {
    if (pow_omega_x_[i] > max_pow_omega_x_) max_pow_omega_x_ = pow_omega_x_[i];
  }

  max_pow_omega_y_ = pow_omega_y_[0];
  for (std::size_t i = 1; i < pow_omega_y_.size(); i++) {
    if (pow_omega_y_[i] > max_pow_omega_y_) max_pow_omega_y_ = pow_omega_y_[i];
  }

  max_pow_x0_ = pow_x0_[0];
  for (std::size_t i = 1; i < pow_x0_.size(); i++) {
    if (pow_x0_[i] > max_pow_x0_) max_pow_x0_ = pow_x0_[i];
  }

  max_pow_y0_ = pow_y0_[0];
  for (std::size_t i = 1; i < pow_y0_.size(); i++) {
    if (pow_y0_[i] > max_pow_y0_) max_pow_y0_ = pow_y0_[i];
  }
}