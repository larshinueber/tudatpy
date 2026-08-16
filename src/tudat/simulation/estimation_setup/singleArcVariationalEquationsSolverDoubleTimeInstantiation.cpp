#include "tudat/simulation/estimation_setup/singleArcVariationalEquationsSolver.h"

namespace tudat
{

namespace propagators
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template class SingleArcVariationalEquationsSolver< double, Time >;
#endif

}  // namespace propagators

}  // namespace tudat
