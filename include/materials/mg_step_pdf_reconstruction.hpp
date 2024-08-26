#ifndef MG_STEP_PDF_RECONSTRUCTION_H
#define MG_STEP_PDF_RECONSTRUCTION_H

#include <materials/mg_angle_distribution.hpp>
#include <utils/error.hpp>

//===========================
// Step-Function PDF 
// P(x) = { 0 for -1 <= x < c and a for c<= x <= 1 } when mean of cosine is positive
// P(x) = { a for -1 <= x < c and 0 for c<= x <= 1 } when mean of cosine is negative
class MGStepReconstruction : public MGAngleDistribution{
    public:
    
    MGStepReconstruction(double p1_moment, std::size_t mat_id);

    double pdf(double x) const override final;

    double cdf(double x) const;

    std::pair<double, double> sample_mu(RNG& rng) const override final;

    
    private:
        double mu_mean_, start_step_point_, step_value_, inv_step_value_;
        bool check_positive_mean = true;

};


#endif