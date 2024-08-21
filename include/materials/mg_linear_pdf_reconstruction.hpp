#ifndef MG_LINEAR_PDF_RECONSTRUCTION_H
#define MG_LINEAR_PDF_RECONSTRUCTION_H

#include <materials/mg_angle_distribution.hpp>
#include <utils/error.hpp>

//=================================
// P(x) =  0 for -1 <= x <= x_0
// P(x) = a + b * x for x_0 <= x <=1
class MGLinearReconstruction : public MGAngleDistribution {
    public:
    MGLinearReconstruction() {
        fatal_error("a linear approximate distribution can be constructed, a cosine mean is required.");
    }

    MGLinearReconstruction(double p1_moment);

    double pdf(double x) const override final;

    double cdf(double x) const;

    std::pair<double, double> sample_mu(RNG& rng) const override final;
    
    private:
        double mu_mean_, interscept_, slope_, x_0, A, B, C_;
        bool check_mu_mean_positive = true;
};

#endif