#include "tudat/simulation/propagation_setup/hybridArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template class HybridArcDynamicsSimulator< double, Time >;
#endif

}  // namespace propagators

}  // namespace tudat
