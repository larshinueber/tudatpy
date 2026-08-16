#include "tudat/simulation/propagation_setup/dynamicsSimulatorBase.h"

namespace tudat
{

namespace propagators
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template class DynamicsSimulator< double, Time >;
#endif

}  // namespace propagators

}  // namespace tudat
