#ifndef TEST_H
#define TEST_H

#include <tallies/zernike_polynomial.hpp>
#include <utils/polynomial_expression.hpp>
#include <utils/error.hpp>
#include <utils/constants.hpp>
#include <utils/zernike_cartesian_polynomial.hpp>
#include <utils/gauss_quadrature.hpp>

#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>

inline void test(){

    std::cout << "========================================" << std::endl;

    Position rp(0.1, 0.44, 0.2);
    Direction u(0.33, 0.22, 0.11);
    const double dist = 0.001;
    Position r = rp + dist * u; 
    

    const double scaled_r = std::sqrt(r.x()*r.x() + r.y()*r.y());
    std::cout << "Checke the radius = " << scaled_r << std::endl;

    double theta = std::atan2(r.y(), r.x());
    if (theta < 0.) {
     theta += 2 * PI;
    }

    ZernikePolynomials zr_poly(11);
    const double value1 = zr_poly.evaluate_zernikes(scaled_r, theta)[11];
    std::cout << std::setprecision(15) << value1 << std::endl;

    ZernikeCartesianPolynomial zr_cart_poly(11);

    std::cout << std::setprecision(15) << zr_cart_poly.evaluate_zernikes(r.x(), r.y())[11] << std::endl;
    std::cout << std::setprecision(15) << zr_cart_poly.evaluate_zernikes(dist, u.x(), u.y(), rp.x(), rp.y())[11] << std::endl;


    std::cout << "\nLet's integrate" << std::endl;
    GaussQuadrature gauss_quad = GaussQuadrature(GaussLegendreQuad<6>());

    std::cout << "size = " << gauss_quad.size() << std::endl;
    double value = 0.;
    for (std::size_t i = 0; i < gauss_quad.size(); i++){
        const auto& absc_and_weight = gauss_quad.quadrature_set()[i];
        const double abscia_quad = absc_and_weight.abscissa;
        const double weight_quad = absc_and_weight.weight;
        const double dist_n_quad = (abscia_quad + 1.) * 0.5 * dist;

        Position r_temp = rp + dist_n_quad * u;
        double scale_r = std::sqrt(r_temp.x()*r_temp.x() + r_temp.y()*r_temp.y());
        theta = std::atan2(r_temp.y(), r_temp.x());
        if (theta < 0.) {
            theta += 2 * PI;
        }
        const double value2 = zr_poly.evaluate_zernikes(scale_r, theta)[11];
        value += 0.5 * dist * value2 * weight_quad;
    }
    std::cout << std::setprecision(15) << "quadrature = " << value << std::endl;
    std::cout << std::setprecision(15) << "Exact int. = ";
    std::cout << zr_cart_poly.line_integrate_zernike(0., dist, u.x(), u.y(), rp.x(), rp.y())[11] << std::endl;

    std::cout << "\n========================================" << std::endl;

    fatal_error("");
}

#endif