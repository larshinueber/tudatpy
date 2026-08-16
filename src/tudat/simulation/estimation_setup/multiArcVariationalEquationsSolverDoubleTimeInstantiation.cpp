#include "tudat/simulation/estimation_setup/multiArcVariationalEquationsSolver.h"

namespace tudat
{

namespace propagators
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template class MultiArcVariationalEquationsSolver< double, Time >;
#endif

}  // namespace propagators

}  // namespace tudat
