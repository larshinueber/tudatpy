#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"

namespace tudat
{

namespace simulation_setup
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template std::shared_ptr< propagators::VariationalEquationsSolver< double, Time > > createVariationalEquationsSolver< double, Time >(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
        const bool integrateEquationsOnCreation );
#endif

}  // namespace simulation_setup

}  // namespace tudat
