#include <materials/mg_step_pdf_reconstruction.hpp>

MGStepReconstruction::MGStepReconstruction(double p1_moment)
    : mu_mean_(p1_moment), start_step_point_(), step_value_()
    {
    if ( mu_mean_ >= 0.0 ){
        check_positive_mean = true;
        }
    else 
        check_positive_mean = false;

    if ( (check_positive_mean== true) && (mu_mean_ == 1.0)){
        std::cout<<"a step approximate cannot be reconstructed, as the first moment is 1.";
        abort();
    }
    if ((check_positive_mean== false) && (mu_mean_ == -1.0)){
        std::cout<<"a step approximate cannot be reconstructed, as the first moment is 1.";
        abort();
    }
    if (check_positive_mean == true){
        start_step_point_ = 2.0 * mu_mean_ - 1.0;
        step_value_ = 0.5 / ( 1.0 - mu_mean_ );
        inv_step_value_ = ( 1.0 - mu_mean_ ) * 2;
    
    }else{
        start_step_point_ = 2.0 * mu_mean_ + 1.0;
        step_value_ = 0.5 / ( 1.0 + mu_mean_ );
        inv_step_value_ = 2 * ( 1.0 + mu_mean_ );
    }
}

double MGStepReconstruction::pdf(double x) const {
    if (check_positive_mean == true){
        if ( x >= start_step_point_ && x <= 1.0 )
            return step_value_;
    
    }else{
        if ( x <= start_step_point_ )
            return step_value_;
    }

    return 0;
}

double MGStepReconstruction::cdf(double x) const {
    if (check_positive_mean==true){
        if ( x >= start_step_point_ && x <= 1.0 )
        return step_value_ * ( x - start_step_point_ );
    
    }else{
        if ( x <= start_step_point_ )
            return step_value_ * ( x + 1.0 );
        else
            return 1.;
    }

    return 0;
}

std::pair<double, double> MGStepReconstruction::sample_mu(RNG& rng) const{ 
    const double random_number = rng();
    // Using the inverse method
    if (check_positive_mean==true){
        if ( random_number == 0.0 ){
            const double scaled_number = rng() * ( start_step_point_ + 1 ) - 1 ;
            return {scaled_number, 1.};
        }
        return {inv_step_value_ * random_number + start_step_point_, 1.};
        
    }else{
        if ( random_number == 1.0){
            const double scaled_number = rng() * ( 1.0 - start_step_point_ ) + start_step_point_;
            return {scaled_number, 1.};
        }
        return {inv_step_value_ * random_number - 1.0, 1.};
    }

    // it should never reach here.
    fatal_error("In step-pdf-reconstruction sample_mu, something weird happened.");
    return {0., 1.};
}