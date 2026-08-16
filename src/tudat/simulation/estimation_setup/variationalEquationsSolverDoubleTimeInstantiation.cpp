#include "tudat/simulation/estimation_setup/variationalEquationsSolverBase.h"

namespace tudat
{

namespace propagators
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS_DOUBLE_TIME
template class VariationalEquationsSolver< double, Time >;
#endif

}  // namespace propagators

}  // namespace tudat
