#include <materials/mg_discrete_angle_distribution.hpp>
#include <boost/math/special_functions/legendre.hpp>

#include <algorithm>
#include <sstream>

MGDiscreteAngleDistribution::MGDiscreteAngleDistribution(const std::vector<double>& legendre_moments, std::size_t mat_id)
    : mu_(), weights_(), cdf_(){
    
    // heighest order 
    const std::size_t L = legendre_moments.size()-1;
    // number of discrete angles
    std::size_t N = (L+1)/2 + 1;
    if ( L % 2 == 0 ){
        N = L/2 + 1;
    }

    // get the roots of N-order legendre polynomial, the returning roots will be only positive,
    // therefore create the negative counter part 
    std::vector<double> roots = boost::math::legendre_p_zeros<double>(static_cast<int>(N));

    // the root of the polynomial will be the discrete angles
    mu_.reserve(N);
    weights_.reserve(N);
    cdf_.reserve(N);
    // add the negative counter-part of angles
    mu_.insert(mu_.end(), roots.rbegin(), roots.rend());
    for ( std::size_t i = 0; i < roots.size(); i++ ){
        mu_[i] *= -1.;    
    }
    if ( N % 2 == 1 ){
        mu_.insert(mu_.end(), roots.begin()+1, roots.end());
    } else {
        mu_.insert(mu_.end(), roots.begin(), roots.end());
    }

    // evaluate the weights
    double normalization_const = 0.;
    for ( std::size_t i = 0; i < mu_.size(); i++){
        const double angle = mu_[i];
        double numerator = 0.;
        for ( std::size_t l = 0; l <= L; l++ ){
            numerator += static_cast<double>(2*l + 1) * 0.5 * legendre_moments[l] * boost::math::legendre_p<double>(static_cast<int>(l), angle);
        }

        double denomenator = 0.;
        for ( std::size_t j = 0; j <= N; j++){
            const double P_j_angle = boost::math::legendre_p<double>(static_cast<int>(j), angle);
            denomenator += static_cast<double>(2*j + 1) * 0.5 * P_j_angle * P_j_angle;
        }
        const double w_i = numerator / denomenator;
        
        if (w_i < 0.){
            std::stringstream mssg;
            mssg << "The weights in the discrete-angle is negative with material id:" <<  mat_id;
            fatal_error(mssg.str());
        }
        weights_.push_back( w_i );
        normalization_const += w_i;
    }

    // cdf
    double cdf0 = 0.;
    cdf_.push_back(cdf0);
    for ( std::size_t i = 0; i < N; i++){
        weights_[i] /= normalization_const;
        cdf0 += weights_[i] ;
        cdf_.push_back(cdf0);
    }
  
}

double MGDiscreteAngleDistribution::pdf(double mu) const {
    auto mu_it = std::lower_bound(mu_.begin(), mu_.end(), mu);
    std::size_t index = static_cast<std::size_t>(std::distance(mu_.begin(), mu_it));
    return weights_[index];
}


std::pair<double, double> MGDiscreteAngleDistribution::sample_mu(RNG& rng) const{
    const double random_number = rng();
    auto cdf_it = std::lower_bound(cdf_.begin(), cdf_.end(), random_number);
    std::size_t index = static_cast<std::size_t>(std::distance(cdf_.begin(), cdf_it));
    return {mu_[index], 1.};
}

