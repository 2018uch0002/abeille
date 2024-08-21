#ifndef MG_LINAR_DELTA_PDF_RECONSTRUCTION_H
#define MG_LINAR_DELTA_PDF_RECONSTRUCTION_H

#include <materials/mg_angle_distribution.hpp>
#include <utils/error.hpp>

//===============================
// Linear-Dirac Delta PDF
// P(x) = a + b * x + c * delta(x-mu_0)
// where mu_0 = 1 if slope is positive
// mu_0 = -1 if slope is negative
class MGLinearDeltaReconstruction : public MGAngleDistribution {
    public:
        MGLinearDeltaReconstruction(){
            fatal_error("a linear-delta approximate distribution can be constructed, a cosine mean is required.");
        }
        
        MGLinearDeltaReconstruction(double p1_moment);

    double pdf(double x) const override final;

    double cdf( double x) const;

    std::pair<double, double> sample_mu(RNG& rng) const override final;

    private:
        // Note that interscept  and slope will be same or opposite in sign.
        double mu_mean_, interscept_, slope_, coeff_delta, mu0, sqrt_2_a; // variable sqrt_2_a = sqrt(2/a)
        bool check_mu_mean_positive = true;

};

#endif