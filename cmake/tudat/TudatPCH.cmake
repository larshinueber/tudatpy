# Shared precompiled-header header list, reused by TUDAT_ADD_LIBRARY(... WITH_PCH)
# across the core tudat libraries. Kept as a single canonical list to avoid
# per-library drift. Scoped to standard-library and Eigen headers only.
set(TUDAT_STD_EIGEN_PCH_HEADERS
        <algorithm>
        <cmath>
        <fstream>
        <functional>
        <iomanip>
        <iostream>
        <limits>
        <map>
        <memory>
        <sstream>
        <string>
        <vector>
        <Eigen/Core>
        <Eigen/Geometry>
        )

# Phase 2: tudat's own environment headers. Scoped to the "setup" libraries that
# actually build/orchestrate Body objects (estimation_setup, environment_setup,
# propagation_setup) rather than all TUDAT_ADD_LIBRARY targets -- body.h and
# createBodiesFactory.h churn more often than the std/Eigen list above, and a
# transitive-include check showed the other core libraries barely reach them
# (0-10%), so PCH-ing this pair there would add rebuild cost for no benefit.
set(TUDAT_ENVIRONMENT_PCH_HEADERS
        <tudat/simulation/environment_setup/body.h>
        <tudat/simulation/environment_setup/createBodiesFactory.h>
        )
