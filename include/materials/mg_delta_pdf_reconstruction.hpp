#ifndef MG_DISCRETE_PDF_RECONSTRUCTION_H
#define MG_DISCRETE_PDF_RECONSTRUCTION_H

#include <materials/mg_angle_distribution.hpp>
#include <utils/error.hpp>
#include <utils/constants.hpp>

//==============================
// Discrete Reconstruction
// P(x) = delta(x - mean) 
class MGDeltaReconstruction : public MGAngleDistribution{
    public:
    MGDeltaReconstruction( double p1_moment) : mu_mean_(p1_moment){}

        double pdf(double /* x */) const override final {
            return INF;
        }

        double cdf(double /* x */ ) {
            return 1.;
        }

        std::pair<double, double> sample_mu(RNG& /*rng*/) const  override final{
            return {mu_mean_, 1.};
        }
    
    private:
    double mu_mean_;

};


#endif