#include "tudat/simulation/estimation_setup/createObservationManager.h"

namespace tudat
{

namespace observation_models
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template std::shared_ptr< ObservationManagerBase< double, Time > > createObservationManagerBase< double, Time >(
        const ObservableType observableType,
        const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettingsList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< Time > > dependentVariablesInterface );

template std::map< ObservableType, std::shared_ptr< ObservationManagerBase< double, Time > > >
createObservationManagersBase< double, Time >(
        const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > fullParameters,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< Time > > dependentVariablesInterface );
#endif

}  // namespace observation_models

}  // namespace tudat
