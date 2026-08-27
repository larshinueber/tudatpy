/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_observations_simulation_settings_bindings.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{
namespace observations_simulation_settings
{

void expose_observation_simulation_settings_factory_bindings( py::module& m )
{
    m.def( "observation_settings_from_collection",
           &tss::getObservationSimulationSettingsFromObservations< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection" ),
           py::arg( "bodies" ),
           R"doc(No documentation found.)doc" );

    m.def( "change_simulation_settings_observable_types",
           &tom::changeObservableTypesOfObservationSimulationSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_simulation_settings" ),
           py::arg_v( "replacement_observable_types",
                      std::map< tom::ObservableType, tom::ObservableType >{
                              { tom::dsn_n_way_averaged_doppler, tom::n_way_differenced_range },
                              { tom::dsn_one_way_averaged_doppler, tom::one_way_differenced_range } },
                      "..." ),
           R"doc(
Change selected observable types in observation simulation settings.

Parameters
----------
observation_simulation_settings : list[ObservationSimulationSettings]
    Observation simulation settings to update.
replacement_observable_types : dict[ObservableType, ObservableType], default = DSN Doppler replacements
    Mapping from each observable type to its replacement. By default, DSN one-way and n-way averaged Doppler types are replaced by their generic averaged Doppler counterparts.
)doc" );
}

}  // namespace observations_simulation_settings
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
