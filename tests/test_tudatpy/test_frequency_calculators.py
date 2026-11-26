"""Check the Python interface of band-specific transmitting frequency calculators.

The behaviour of these calculators is tested in C++ (unitTestTransmittingFrequenciesMultiBand.cpp);
these tests only check that the interface is exposed and that errors reach Python.
"""

import pytest

from tudatpy.dynamics import environment, environment_setup
from tudatpy.dynamics.environment import FrequencyBands


def constant_ramp(frequency: float) -> environment.PiecewiseLinearFrequencyInterpolator:
    """Ramp table with a single constant-frequency ramp from 0 to 100 s."""
    return environment.PiecewiseLinearFrequencyInterpolator([0.0], [100.0], [0.0], [frequency])


def test_multi_band_frequency_interpolator():
    """A multi-band calculator keeps a separate frequency calculator for each band."""
    calculator = environment.MultiBandFrequencyInterpolator(
        {FrequencyBands.x_band: constant_ramp(7.18e9)}
    )
    calculator.add_band_interpolator(FrequencyBands.ka_band, constant_ramp(34.3e9))

    assert calculator.has_band(FrequencyBands.ka_band)
    assert not calculator.has_band(FrequencyBands.s_band)
    assert sorted(calculator.available_bands, key=int) == [
        FrequencyBands.x_band,
        FrequencyBands.ka_band,
    ]
    assert list(calculator.get_band_interpolator(FrequencyBands.ka_band).start_frequencies) == [
        34.3e9
    ]


def test_band_specific_transmitting_frequency_calculators():
    """Ground stations and vehicles keep a separate transmitting frequency calculator for each band."""
    settings = environment_setup.BodyListSettings("SSB", "J2000")
    settings.add_empty_settings("Earth")
    settings.get("Earth").shape_settings = environment_setup.shape.spherical(6378137.0)
    bodies = environment_setup.create_system_of_bodies(settings)
    environment_setup.add_ground_station(
        bodies.get("Earth"),
        environment_setup.ground_station.basic_station("DSS-43", [6378137.0, 0.0, 0.0]),
    )
    station = bodies.get("Earth").get_ground_station("DSS-43")

    station.add_transmitting_frequency_calculator(constant_ramp(7.18e9), FrequencyBands.x_band)
    station.add_transmitting_frequency_calculator(constant_ramp(2.11e9), FrequencyBands.s_band)
    assert isinstance(
        station.transmitting_frequency_calculator, environment.MultiBandFrequencyInterpolator
    )
    assert station.has_frequency_calculator(FrequencyBands.s_band)
    assert not station.has_frequency_calculator(FrequencyBands.ka_band)
    s_band_calculator = station.get_transmitting_frequency_calculator(FrequencyBands.s_band)
    assert list(s_band_calculator.start_frequencies) == [2.11e9]

    # Band-specific and band-independent calculators cannot be combined
    with pytest.raises(RuntimeError):
        station.add_transmitting_frequency_calculator(
            environment.ConstantTransmittingFrequencyCalculator(7.18e9)
        )

    vehicle_systems = environment.VehicleSystems()
    vehicle_systems.add_transmitted_frequency_calculator(
        constant_ramp(7.18e9), FrequencyBands.x_band
    )

    sband_to_add = constant_ramp(2.11e9)
    with pytest.raises(RuntimeError):
        vehicle_systems.add_transmitted_frequency_calculator(sband_to_add)
