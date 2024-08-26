#include <materials/mg_linear_delta_pdf_reconstruction.hpp>
#include <utils/constants.hpp>

MGLinearDeltaReconstruction::MGLinearDeltaReconstruction(double p1_moment)
    : mu_mean_(p1_moment){

    if ( mu_mean_ >= 0.0 ){
        check_mu_mean_positive = true;
        interscept_ = 3.0 * 0.25 * (1.0 - mu_mean_);
        slope_ = interscept_;
        coeff_delta = 0.5 * ( 3.0 * mu_mean_ - 1.0);
        mu0 = 1.0;
    }

    if (mu_mean_ < 0.0){
        check_mu_mean_positive = false;
        interscept_ = 3.0 * 0.25 * (1 + mu_mean_);
        slope_ = -1.0 * interscept_;
        coeff_delta = 1.0 - 2 * interscept_;
        mu0 = -1.0;
    }

    sqrt_2_a = std::sqrt( 2/interscept_ );

}

double MGLinearDeltaReconstruction::pdf(double x) const {
    double pdf_value = interscept_ + slope_ * x;

    if ( x == mu0 ){
        return INF;
    }
    return pdf_value;
}


double MGLinearDeltaReconstruction::cdf( double x) const {
    double cdf_value = interscept_ * (x+1) + slope_ * 0.5 * (x*x - 1.0);

    if ( mu0 == 1.0){
        if ( x == mu0 )
            return 1.0;
        return cdf_value;
    }

    if ( mu0 == -1.0){
        if ( x == mu0 )
            return 0.0;

        return cdf_value + coeff_delta;
    }
    return 0.;
}

std::pair<double, double> MGLinearDeltaReconstruction::sample_mu(RNG& rng) const {
    const double random_number = rng();
    if (mu0 == 1.0){
        if( random_number > (1-coeff_delta)) 
            return {1., 1.};
        
        const double sampled_mu = -1 + sqrt_2_a * std::sqrt(random_number);
        return {sampled_mu, 1.};
    }

    if (mu0 == -1.0){
        if ( random_number <= coeff_delta)
            return {-1., 1.};
        const double sampled_mu = 1 - sqrt_2_a * std::sqrt(1-random_number);
        return {sampled_mu, 1.};
    }

     // it should never reach here.
    fatal_error("In linear-delta-pdf-reconstruction sample_mu, something weird happened.");
    return {mu0, 1.};

}