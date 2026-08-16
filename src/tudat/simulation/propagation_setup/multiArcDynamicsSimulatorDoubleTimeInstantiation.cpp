#include "tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template class MultiArcDynamicsSimulator< double, Time >;
#endif

}  // namespace propagators

}  // namespace tudat
