import sys

from tudatpy._deprecation import register_deprecated_object
from tudatpy.kernel.estimation.observations_setup.observations_simulation_settings import *
from tudatpy.kernel.estimation.observable_models.observables_simulation import (
    create_observation_simulators,
)

register_deprecated_object(
    sys.modules[__name__],
    "create_observation_simulators",
    create_observation_simulators,
    new_name="tudatpy.estimation.observable_models.observables_simulation.create_observation_simulators",
)
