#include <materials/mg_semicontinuous_distribution.hpp>
#include <utils/error.hpp>
#include <utils/constants.hpp>

#include <cmath>
#include <sstream>


MGSemicontinuousDistribution::MGSemicontinuousDistribution(const std::vector<double> legendre_moments, std::size_t mat_id)
    : x1_(), x2_(), p_(), beta_(), M1_(), M2_(), M3_(), p1_pdf_negative(false){

    // check the legendre-moments of three orders is there or not
    if ( legendre_moments.size() < 4 ){
        std::stringstream mssg;
        mssg<<"The moments must be at least three for semicontinuous scatter for material id: "<< mat_id ;
        fatal_error(mssg.str());
    }

    M1_ = legendre_moments[1];
    M2_ = (1. + 2. * legendre_moments[2])/3.;
    M3_ = (3. * legendre_moments[1] + 2. * legendre_moments[3]) / 5.;

    M1_star_ = M1_;
    M2_star_ = 1./3. ;
    M3_star_ = 0.6 * M1_;

    // the semi-continuous has a line in the background, which is constructed from the p1-moments, therefore, if the pdf consutructed using the p1 moments, may come as a negative. To avoid this a linear + delta at(mu = 1) while conserving the first-moment is used. 
    if ( std::abs(M1_*3.) > 1. ){
        p1_pdf_negative = true;
        M2_star_ = M1_;
        M3_star_ = 0.2 * (6.*M1_ - 1.);
    }

    // evaluate all the intermediate variable to get the beta_
    const double sigma_1_square = M2_ - M1_*M1_;
    const double L2 = M3_ - M2_*M1_;

    const double sigma_1_square_star = M2_star_ - M1_star_ * M1_star_;
    const double L2_star = M3_star_ - M2_star_ * M1_star_;

    beta_ = 1.;

    const double A_plus = (1. - M2_) * sigma_1_square / ( 1 - M1_ ) - L2;
    const double a_plus = sigma_1_square / std::sqrt(1. - M1_);
    const double A_plus_star = (1. - M2_star_) * sigma_1_square_star / ( 1 - M1_star_ ) - L2_star;
    const double a_plus_star = sigma_1_square_star / std::sqrt(1. - M1_star_);
    const double B1_plus = A_plus + A_plus_star + (a_plus - a_plus_star)*(a_plus - a_plus_star);
    const double term1_plus = (A_plus - A_plus_star)*(A_plus - A_plus_star);
    const double term2_plus = (a_plus - a_plus_star)*(a_plus - a_plus_star)*(a_plus - a_plus_star)*(a_plus - a_plus_star); 
    const double term3_plus = 2. * (A_plus + A_plus_star) * (a_plus - a_plus_star)*(a_plus - a_plus_star);
    const double B2_plus = std::sqrt( term1_plus + term2_plus + term3_plus );
    const double beta1_plus = (B1_plus + B2_plus) * 0.5 / A_plus_star;
    if ( beta1_plus >= 0. && beta1_plus <= 1.){
        if ( beta_ > beta1_plus )
            beta_ = beta1_plus;
    }
    const double beta2_plus = (B1_plus - B2_plus) * 0.5 / A_plus_star;
    if ( beta2_plus >= 0. && beta2_plus <= 1. ){
        if ( beta_ > beta2_plus )
            beta_ = beta2_plus;
    }

    const double A_minus = (1. - M2_) * sigma_1_square / ( 1 + M1_ ) + L2;
    const double a_minus = sigma_1_square / std::sqrt(1. + M1_);
    const double A_minus_star = (1. - M2_star_) * sigma_1_square_star / ( 1 + M1_star_ ) + L2_star;
    const double a_minus_star = sigma_1_square_star / std::sqrt(1. + M1_star_);
    const double B1_minus = A_minus + A_minus_star + (a_minus - a_minus_star)*(a_minus - a_minus_star);
    const double term1_minus = (A_minus - A_minus_star)*(A_minus - A_minus_star);
    const double term2_minus = (a_minus - a_minus_star)*(a_minus - a_minus_star)*(a_minus - a_minus_star)*(a_minus - a_minus_star); 
    const double term3_minus = 2. * (A_minus + A_minus_star) * (a_minus - a_minus_star)*(a_minus - a_minus_star);
    const double B2_minus = std::sqrt( term1_minus + term2_minus + term3_minus );
    const double beta1_minus = (B1_minus + B2_minus) * 0.5 / A_minus_star;
    if ( beta1_minus >= 0. && beta1_minus <= 1.){
        if ( beta_ > beta1_minus)
            beta_ = beta1_minus;
    }
    const double beta2_minus = (B1_minus - B2_minus) * 0.5 / A_minus_star;
    if ( beta2_minus >= 0. && beta2_minus <= 1. ){
        if ( beta_ > beta2_minus)
            beta_ = beta2_minus;
    }


    const double L2_bar = (L2 - beta_ * L2_star) / (1.-beta_);
    const double sigma_1_sq_bar = (sigma_1_square - beta_ * sigma_1_square_star) / (1. - beta_);

    const double lambda_bar = L2_bar / sigma_1_sq_bar;
    const double S_bar = std::sqrt( (lambda_bar - 2. * M1_)*(lambda_bar - 2. * M1_) + 4.*sigma_1_sq_bar );

    x1_ = 0.5 * (lambda_bar + S_bar);
    x2_ = 0.5 * (lambda_bar - S_bar);

    p_ = 2. * sigma_1_sq_bar / ( S_bar*S_bar + S_bar*(lambda_bar - 2.*M1_) );

}

double MGSemicontinuousDistribution::pdf(double x) const {
    // if pdf constructed by p1 is positive
    if ( p1_pdf_negative == false){
        const double linear_part = 0.5 * (1 + 3 * M1_ * x) * beta_;
        if ( x == x1_ || x == x2_ ){
            return INF;
        }
        return linear_part;
    }
    // if pdf constructed by p1 is negative
    const double linear_part = 0.25 * 3. * (1.- M1_) * (1. + x);
    if ( x == 1.){
        return INF;
    }
    return linear_part;
}

std::pair<double, double> MGSemicontinuousDistribution::sample_mu(RNG& rng) const {
    const double xi = rng();
    if ( xi < beta_ ){
        double a = 0.5;
        double b = 1.5 * M1_;
        const double xi_2 = rng();
        
        // if pdf constructed by p1 is positive
        if ( p1_pdf_negative == false){
            return {(-a + std::sqrt((a-b)*(a-b) + 2. * b * xi_2)) / b , 1.};
        } else {
            if ( xi_2 <= 0.5 * (3*M1_ - 1.)){
                return {1., 1.};
            } 
            return {-1. + std::sqrt(2 * xi_2 / (0.75 * (1.-M1_))), 1.};
        }
    } else {
        const double xi_3 = rng();
        if ( xi_3 < p_){
            return {x1_, 1.};
        } else {
            return {x2_, 1.};
        }
    }

    return {xi, 1.};
}