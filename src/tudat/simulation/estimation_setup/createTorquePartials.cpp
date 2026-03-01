#include "tudat/astro/orbit_determination/rotational_dynamics_partials/constantTorquePartial.h"
#include "tudat/simulation/estimation_setup/createTorquePartials.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create torque partial to be used for constant torques in angular acceleration
std::shared_ptr< acceleration_partials::TorquePartial > createConstantTorqueRotationalDynamicsPartial(
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const basic_astrodynamics::SingleBodyTorqueModelMap& torqueVector )
{
    std::function< Eigen::Vector3d( ) > angularVelocityFunction =
            [body = acceleratedBody.second]( ) { return body->getCurrentAngularVelocityVectorInLocalFrame( ); };
    std::function< Eigen::Matrix3d( ) > inertiaTensorFunction =
            [body = acceleratedBody.second]( ) { return body->getBodyInertiaTensor( ); };

    std::function< double( ) > inertiaTensorNormalizationFunction;
    if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( acceleratedBody.second->getGravityFieldModel( ) ) !=
        nullptr )
    {
        inertiaTensorNormalizationFunction =
                [shGravityField = std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                         acceleratedBody.second->getGravityFieldModel( ) )]( ) { return shGravityField->getInertiaTensorNormalizationFactor( ); };
    }

    std::function< double( ) > gravitationalParameterFunction;
    if( acceleratedBody.second->getGravityFieldModel( ) != nullptr )
    {
        gravitationalParameterFunction =
                [gravityModel = acceleratedBody.second->getGravityFieldModel( )]( ) { return gravityModel->getGravitationalParameter( ); };
    }

    return std::make_shared< acceleration_partials::ConstantTorquePartial >( angularVelocityFunction,
                                                                             inertiaTensorFunction,
                                                                             inertiaTensorNormalizationFunction,
                                                                             gravitationalParameterFunction,
                                                                             torqueVector,
                                                                             acceleratedBody.first );
}

}  // namespace simulation_setup

}  // namespace tudat
