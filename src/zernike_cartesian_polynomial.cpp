#include <utils/zernike_cartesian_polynomial.hpp>

ZernikeCartesianPolynomial::ZernikeCartesianPolynomial(const std::size_t order)
    : zr_poly_expression_(), order_(order), max_n_() {
  // Evaluate the polynomial coefficients up to given order
  // zr_poly_expression_ will have the polynomial-expression stating from lowest
  // power to higher for each order.
  zr_poly_expression_.reserve(order_ + 1);

  // get a maximum size to store the coefficients
  std::pair<std::size_t, int> high_n_and_l = get_n_and_l(order_);
  // store the n for highest order
  max_n_ = high_n_and_l.first;
  const int l_max_n = high_n_and_l.second;
  const std::size_t m_max_n = static_cast<std::size_t>(std::abs(l_max_n));

  // evaluate the polynomial exression of cosine and sine for upto m*theta
  const PolynomialExpression r_cosine("X");
  const PolynomialExpression r_sine("Y");
  std::vector<PolynomialExpression> r_cosine_m_theta;
  r_cosine_m_theta.reserve(max_n_);
  r_cosine_m_theta.push_back(r_cosine);

  std::vector<PolynomialExpression> r_sine_m_theta;
  r_sine_m_theta.reserve(max_n_);
  r_sine_m_theta.push_back(r_sine);

  for (std::size_t it = 1; it < max_n_; it++) {
    // cos (m * theta) = cos( (m-1)*theta ) * cos( theta ) - sin( (m-1)*theta )
    // * sin( theta )
    const PolynomialExpression cos_m_theta =
        r_cosine_m_theta[it - 1] * r_cosine - r_sine_m_theta[it - 1] * r_sine;
    r_cosine_m_theta.push_back(cos_m_theta);

    // sin (m * theta) = sin (m-1)*theta ) * cos( theta ) + cos( (m-1)*theta ) *
    // sin( theta )
    const PolynomialExpression sin_m_theta =
        r_sine_m_theta[it - 1] * r_cosine + r_cosine_m_theta[it - 1] * r_sine;
    r_sine_m_theta.push_back(sin_m_theta);
  }

  // evaluate the polynomial expression for the powers of r_square
  const PolynomialExpression r_square = r_cosine * r_cosine + r_sine * r_sine;
  std::vector<PolynomialExpression> r_square_power;
  const std::size_t size_r_sq_pow =
      static_cast<std::size_t>((max_n_ - m_max_n) / 2) + 1;
  r_square_power.reserve(size_r_sq_pow);
  r_square_power.emplace_back(1., 0);
  r_square_power.push_back(r_square);
  PolynomialExpression r_sqr_k = r_square;
  for (std::size_t it = 1; it < size_r_sq_pow; it++) {
    r_sqr_k *= r_square;
    r_square_power.push_back(r_sqr_k);
  }

  // evaluate the radial polynomial, then multiply it by cosine or sine to
  // construce the Zernike polynomials in cartesian coordinates
  int l;
  std::size_t n, m;
  zr_poly_expression_.emplace_back(1., 0);  // store the 0-th order

  for (std::size_t it = 1; it <= order_; it++) {
    std::pair<std::size_t, int> n_and_l = get_n_and_l(it);
    n = n_and_l.first;
    l = n_and_l.second;
    m = static_cast<std::size_t>(std::abs(l));
    std::size_t max_k = static_cast<std::size_t>((n - m) / 2);

    double sign_change = 1.;
    if (static_cast<int>(max_k % 2) == 1) sign_change = -1.;

    std::size_t k;
    PolynomialExpression poly_ex(0., 0);
    for (std::size_t i = 0; i <= max_k; i++) {
      k = max_k - i;
      const double numeriator = factorial(n - k);
      const std::size_t avg_n_m = (n + m) / 2;
      const std::size_t mid_n_m = (n - m) / 2;
      const double denominator =
          factorial(k) * factorial(avg_n_m - k) * factorial(mid_n_m - k);

      const double coeff = sign_change * numeriator / denominator;
      poly_ex += coeff * r_square_power[i];
      //   if (it == 4){
      //     std::cout << "k = " << k << "\tsign_change" << sign_change <<
      //     "\tcoff = " << coeff << std::endl;
      //   }

      sign_change = -1 * sign_change;  // for the next iteration
    }

    if (l < 0) {
      // polynmial is odd, therefore use the sin( m * theta)
      poly_ex *= r_sine_m_theta[m - 1];
    } else if (l > 0) {
      // polynmial is even, therefore use the cos( m * theta)
      poly_ex *= r_cosine_m_theta[m - 1];
    }

    zr_poly_expression_.push_back(poly_ex);
  }
}

std::vector<double> ZernikeCartesianPolynomial::evaluate_zernikes(
    const double t, const double omega_x, const double omega_y, const double x0,
    const double y0) const {
  std::vector<double> values;
  values.reserve(order_ + 1);
  for (std::size_t it = 0; it <= order_; it++) {
    values.push_back(
        zr_poly_expression_[it].evaluate(t, omega_x, omega_y, x0, y0));
  }
  return values;
}

std::vector<double> ZernikeCartesianPolynomial::evaluate_zernikes(
    const double x0, const double y0) const {
  std::vector<double> values;
  values.reserve(order_ + 1);
  for (std::size_t it = 0; it <= order_; it++) {
    values.push_back(zr_poly_expression_[it].evaluate_at_point(x0, y0));
  }
  return values;
}

std::vector<double> ZernikeCartesianPolynomial::line_integrate_zernike(
    const double t_a, const double t_b, const double omega_x,
    const double omega_y, const double x0, const double y0) const {
  std::vector<double> values;
  values.reserve(order_ + 1);
  for (std::size_t it = 0; it <= order_; it++) {
    values.push_back(
        zr_poly_expression_[it].integrate(t_a, t_b, omega_x, omega_y, x0, y0));
  }
  return values;
}