/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/io/trackingSupplementaryData.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/createObservationCollection.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ground_stations;

BOOST_AUTO_TEST_SUITE( test_frequency_supplementary_data_in_bodies )

using FrequencyDataMap =
        std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::FrequencySupplementaryData > > >;
using FrequencyRamp = data::RampedFrequencySupplementaryData::FrequencyRamp;

//! Create bodies with either a ground station (Earth/DSS-43) or a spacecraft as link end, and return the link end.
std::pair< std::string, std::string > createLinkEndBody( simulation_setup::SystemOfBodies& bodies, const bool useGroundStation )
{
    const std::string bodyName = useGroundStation ? "Earth" : "Spacecraft";
    const std::string referencePointName = useGroundStation ? "DSS-43" : "";
    bodies.createEmptyBody( bodyName );
    if( useGroundStation )
    {
        simulation_setup::createGroundStation( bodies.at( bodyName ), referencePointName, Eigen::Vector3d( 6378137.0, 0.0, 0.0 ) );
    }
    return { bodyName, referencePointName };
}

//! Transmitting frequency calculator of the ground station or spacecraft created with createLinkEndBody.
std::shared_ptr< StationFrequencyInterpolator > getTransmittingCalculator( const simulation_setup::SystemOfBodies& bodies,
                                                                           const bool useGroundStation )
{
    return useGroundStation ? bodies.at( "Earth" )->getGroundStation( "DSS-43" )->getTransmittingFrequencyCalculator( )
                            : bodies.at( "Spacecraft" )->getVehicleSystems( )->getTransmittedFrequencyCalculator( );
}

// Install a band-independent frequency schedule on both a ground station and a spacecraft, then
// verify that each requested frequency starts at the right time and stays constant.
BOOST_AUTO_TEST_CASE( testPiecewiseConstantFrequencySupplementaryDataUsesZeroRateRamps )
{
    const std::map< double, double > firstHistory = { { 0.0, 100.0 }, { 10.0, 150.0 } };
    const std::map< double, double > secondHistory = { { 10.0, 200.0 }, { 20.0, 300.0 } };

    for( const bool useGroundStation : { true, false } )
    {
        simulation_setup::SystemOfBodies bodies;
        FrequencyDataMap frequencyData;
        frequencyData[ createLinkEndBody( bodies, useGroundStation ) ] = {
            std::make_shared< data::PiecewiseConstantFrequencySupplementaryData >( firstHistory ),
            std::make_shared< data::PiecewiseConstantFrequencySupplementaryData >( secondHistory )
        };
        observation_models::setFrequencySupplementaryDataInBodies( bodies, frequencyData );

        // A spacecraft must receive its transmitter information even when it had no systems yet.
        if( !useGroundStation )
        {
            BOOST_REQUIRE( bodies.at( "Spacecraft" )->getVehicleSystems( ) != nullptr );
        }
        auto calculator =
                std::dynamic_pointer_cast< PiecewiseLinearFrequencyInterpolator >( getTransmittingCalculator( bodies, useGroundStation ) );
        BOOST_REQUIRE( calculator != nullptr );

        // The three frequency steps must be constant; the second input supplies the value at time 10.
        BOOST_CHECK( calculator->getRampRates( ) == std::vector< double >( 3, 0.0 ) );
        BOOST_CHECK_EQUAL( calculator->getTemplatedCurrentFrequency<>( 5.0 ), 100.0 );
        BOOST_CHECK_EQUAL( calculator->getTemplatedCurrentFrequency<>( 15.0 ), 200.0 );
        BOOST_CHECK_EQUAL( calculator->getTemplatedCurrentFrequency<>( 25.0 ), 300.0 );
    }
}

// Install band-specific frequency tables on a ground station and a spacecraft in one call, and verify that the tables
// are grouped per band into one calculator per band, and that tables of the same band are merged.
BOOST_AUTO_TEST_CASE( testFrequencySupplementaryDataPerBand )
{
    using namespace tudat::basic_astrodynamics;

    const double xBandFrequency = 7.18322247E9;
    const double xBandRampRate = -45.0;
    const double sBandFrequency = xBandFrequency * 221.0 / 749.0;

    for( const bool useGroundStation : { true, false } )
    {
        simulation_setup::SystemOfBodies bodies;

        // Two X-band ramp tables covering separate time spans, and an S-band table of piecewise constant frequencies
        FrequencyDataMap frequencyData;
        frequencyData[ createLinkEndBody( bodies, useGroundStation ) ] = {
            std::make_shared< data::RampedFrequencySupplementaryData >(
                    std::vector< FrequencyRamp >{ FrequencyRamp( 0.0, 100.0, xBandFrequency, xBandRampRate ) }, "X-band" ),
            std::make_shared< data::PiecewiseConstantFrequencySupplementaryData >(
                    std::map< double, double >{ { 0.0, sBandFrequency }, { 150.0, sBandFrequency + 1.0E3 } }, "S-band" ),
            std::make_shared< data::RampedFrequencySupplementaryData >(
                    std::vector< FrequencyRamp >{ FrequencyRamp( 200.0, 300.0, xBandFrequency + 1.0E3, xBandRampRate ) }, "X-band" )
        };
        observation_models::setFrequencySupplementaryDataInBodies( bodies, frequencyData );

        // The link end must receive a separate calculator for each band
        auto calculator =
                std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >( getTransmittingCalculator( bodies, useGroundStation ) );
        BOOST_REQUIRE( calculator != nullptr );
        BOOST_CHECK( calculator->getAvailableBands( ) ==
                     ( std::vector< FrequencyBands >{ FrequencyBands::s_band, FrequencyBands::x_band } ) );

        // Both X-band tables must be used by the X-band calculator
        BOOST_CHECK_CLOSE_FRACTION( ( calculator->getTemplatedCurrentFrequency< double, double >( 50.0, FrequencyBands::x_band ) ),
                                    xBandFrequency + 50.0 * xBandRampRate,
                                    1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( ( calculator->getTemplatedCurrentFrequency< double, double >( 250.0, FrequencyBands::x_band ) ),
                                    xBandFrequency + 1.0E3 + 50.0 * xBandRampRate,
                                    1.0E-14 );

        // The S-band frequencies must be constant per step
        BOOST_CHECK_EQUAL( ( calculator->getTemplatedCurrentFrequency< double, double >( 100.0, FrequencyBands::s_band ) ),
                           sBandFrequency );
        BOOST_CHECK_EQUAL( ( calculator->getTemplatedCurrentFrequency< double, double >( 200.0, FrequencyBands::s_band ) ),
                           sBandFrequency + 1.0E3 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
