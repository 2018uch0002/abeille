#ifndef MG_DISCRETE_ANGLE_DISTRIBUTION_H
#define MG_DISCRETE_ANGLE_DISTRIBUTION_H

#include <materials/mg_angle_distribution.hpp>

class MGDiscreteAngleDistribution : public MGAngleDistribution{
    public:
        MGDiscreteAngleDistribution(const std::vector<double>& legendre_moments, std::size_t mat_id);

        double pdf(double mu) const override final;

        std::pair<double, double> sample_mu(RNG& rng) const override final; 

    private:
        std::vector<double> mu_;
        std::vector<double> weights_;
        std::vector<double> cdf_;
};

#endif