#include "tudat/simulation/estimation_setup/hybridArcVariationalEquationsSolver.h"

namespace tudat
{

namespace propagators
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template class HybridArcVariationalEquationsSolver< double, Time >;
#endif

}  // namespace propagators

}  // namespace tudat
