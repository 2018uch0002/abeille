#ifndef MG_SEMICONTINUOUS_DISTRIBUTION_H
#define MG_SEMICONTINUOUS_DISTRIBUTION_H

#include <materials/mg_angle_distribution.hpp>

#include <vector>

class MGSemicontinuousDistribution : public MGAngleDistribution {
    public:
    MGSemicontinuousDistribution(const std::vector<double> legendre_moments, std::size_t mat_id);

    double pdf(double mu) const override final;
    
    std::pair<double, double> sample_mu(RNG& rng) const override final;

    private:
    double x1_, x2_, p_, beta_; 
    double M1_, M2_, M3_;
    double M1_star_, M2_star_, M3_star_;
    bool p1_pdf_negative = false;

};

#endif