#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template class SingleArcDynamicsSimulator< double, Time >;
#endif

}  // namespace propagators

}  // namespace tudat
