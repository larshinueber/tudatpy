#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"

namespace tudat
{

namespace simulation_setup
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template std::shared_ptr< propagators::DynamicsSimulator< double, Time > > createDynamicsSimulator< double, Time >(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated );
#endif

}  // namespace simulation_setup

}  // namespace tudat
