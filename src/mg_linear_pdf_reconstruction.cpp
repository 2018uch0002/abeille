#include <materials/mg_linear_pdf_reconstruction.hpp>

MGLinearReconstruction::MGLinearReconstruction(double p1_moment)
    :mu_mean_(p1_moment), interscept_(), slope_(), x_0(), A(), B() {

    // make the flag variable for the postive and negative mean true or false respectively
    if ( mu_mean_ >= 0.0)
        check_mu_mean_positive = true;
    else
        check_mu_mean_positive = false;

    if (mu_mean_ == 1.0 && check_mu_mean_positive == true ){
        fatal_error("a linear approximate cannot be reconstructed, as the first moment is 1.");
    }
    
    if (mu_mean_ == -1.0 && check_mu_mean_positive == false ){
        fatal_error("a linear approximate cannot be reconstructed, as the first moment is -1.");
    }

    if (check_mu_mean_positive == true){
        slope_ = 2.0 / ( 9.0 * (mu_mean_-1)*(mu_mean_-1) );
        x_0 = ( 3.0 * mu_mean_ - 2.0);
        interscept_ = - slope_ * x_0;

        A = slope_ * 0.5;
        B = interscept_ ;
        C_ = -interscept_ * x_0 - 0.5 * slope_ * x_0*x_0;
    }
    else{
        slope_ = -2.0 / ( 9.0 * (mu_mean_+1) * (mu_mean_+1) );
        x_0 = ( 3.0 * mu_mean_ + 2.0 );
        interscept_ = - slope_ * x_0;

        A = 0.5 * slope_;
        B = interscept_;
        C_ = interscept_ - slope_ * 0.5 ;
    }
}

double MGLinearReconstruction::pdf(double x) const {

    if ( check_mu_mean_positive == true ){
        if ( x >= x_0 && x <= 1.0 )
            return ( interscept_ + slope_ * x );
    }
    else{
        if ( x >= -1.0 && x<= x_0)
            return ( interscept_ + slope_ * x);
    }
    return 0;
}

double MGLinearReconstruction::cdf(double x) const {
    if (check_mu_mean_positive== true) {        
        if ( x >= x_0 && x <= 1.0 )
            return ( interscept_*(x - x_0) + 0.5 * slope_ * (x*x - x_0*x_0) );
    }
    else{
        if ( x >= -1.0 && x <= x_0 )
            return ( interscept_* (x + 1.0 ) + 0.5 * slope_ * ( x*x - 1 ) );
        else
        return 1.0;
    }
    return 0;
}

std::pair<double, double> MGLinearReconstruction::sample_mu(RNG& rng) const {
    const double random_number = rng();
    // Using the inverse method
    if ( check_mu_mean_positive == true ){
        if ( random_number == 0.0 ){
            const double scaled_number = rng() * ( x_0 + 1 ) - 1 ;
            return {scaled_number, 1.};
        }
        const double C = (C_- random_number);
        const double discriminant = B*B - 4 * A * C;
        if (discriminant<0.0)
            fatal_error("In-linear reconstruction, while sampling the discriminant is negative.");
        
        const double root1 = ( - B + std::sqrt(discriminant) ) * 0.5 / A;
        if ( root1 < x_0 && root1 >= 1.0 )
            return sample_mu( rng );
        
        return {root1, 1.};
    }
    else{
        if ( random_number == 1.0 ){
            const double scaled_number = rng() * ( 1 - x_0 ) + x_0;
            return {scaled_number, 1.};
        }
        const double C = ( C_ - random_number);
        const double discriminant =  B*B - 4* A* C;
        if ( discriminant < 0.0)
            fatal_error("In-linear reconstruction, while sampling the discriminant is negative.");

        const double root1 = ( -B + std::sqrt(discriminant) ) * 0.5 / A;
        return {root1, 1.};
    }

    // it should never reach here.
    fatal_error("In linear-pdf-reconstruction sample_mu, something weird happened.");
    return {0., 1.};
}