#include <materials/legendre_distribution.hpp>
#include <materials/mg_tabulated_distribution.hpp>
#include <materials/mg_step_pdf_reconstruction.hpp>
#include <materials/mg_linear_pdf_reconstruction.hpp>
#include <materials/mg_linear_delta_pdf_reconstruction.hpp>
#include <materials/mg_delta_pdf_reconstruction.hpp>
#include <materials/mg_discrete_angle_distribution.hpp>
#include <materials/mg_semicontinuous_distribution.hpp>
#include <utils/settings.hpp>
#include <utils/error.hpp>


std::shared_ptr<MGAngleDistribution> make_mg_angle_distribution(const LegendreDistribution& legendre_dist, const std::size_t mat_id){
    
    // if tabulated is distribution type.
    if ( settings::scatter_distribution_type == settings::MGScatterDistributionType::Tabulated){
        return legendre_dist.linearize();
    }

    // note that the legendre_moments are alredy multiplied by the (2*l+1)/2.
    const std::vector<double> legendre_moments = legendre_dist.a();
    std::vector<double> moments ;
    moments.reserve(legendre_moments.size());
    for (std::size_t i = 0; i < legendre_moments.size(); i++){
        moments.push_back(legendre_moments[i] * 2. / static_cast<double>(2*i + 1));
    }

    if (settings::scatter_distribution_type == settings::MGScatterDistributionType::DiscreteAngle){
        return std::make_shared<MGDiscreteAngleDistribution>(moments, mat_id);

    } if (settings::scatter_distribution_type == settings::MGScatterDistributionType::SemiContinuousLinear){
        return std::make_shared<MGSemicontinuousDistribution>(moments, mat_id);

    }else {
        if ( legendre_moments.size() == 1 ){
            std::stringstream mssg;
            mssg << "The approximate distribution cannot be made in the isotropic distribution with material id " << mat_id;
            fatal_error(mssg.str());
        }
        if ( settings::scatter_distribution_type == settings::MGScatterDistributionType::StepApproximate){
            return std::make_shared<MGStepReconstruction>(moments[1], mat_id);
        } else if ( settings::scatter_distribution_type == settings::MGScatterDistributionType::LinearApproximate){
            return std::make_shared<MGLinearReconstruction>(moments[1], mat_id);
        } else if ( settings::scatter_distribution_type == settings::MGScatterDistributionType::LinearDeltaApproximate){
            return std::make_shared<MGLinearDeltaReconstruction>(moments[1]);
        } else if (settings::scatter_distribution_type == settings::MGScatterDistributionType::DeltaApproximate){
            return std::make_shared<MGDeltaReconstruction>(moments[1]);
        }
    }

    // it shouldn't reach here, a unknow type of scattering-distribution is given
    std::stringstream mssg;
    mssg << "unknown type of scattering distribtion is given with material id " << mat_id;
    fatal_error(mssg.str());
    return legendre_dist.linearize();
}