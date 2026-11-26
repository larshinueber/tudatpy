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
#include <set>
#include <string>
#include <vector>

#include "tudat/io/preProcessOdfFile.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/createObservationCollection.h"
#include "tudat/simulation/estimation_setup/compressDopplerObservationCollection.h"

using namespace tudat::input_output;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;

using namespace tudat;

BOOST_AUTO_TEST_SUITE( test_process_odf_data )

namespace
{

const std::string syntheticDopplerStationName = "DSS-14";

double getSyntheticDopplerStartUtc( )
{
    return basic_astrodynamics::DateTime( 2000, 1, 1, 12, 0, 0.0 ).epoch< double >( );
}

std::vector< double > convertUtcOffsetsToTdbTimes( const std::vector< double >& utcOffsets )
{
    earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter;
    Eigen::Vector3d stationPosition = simulation_setup::getCombinedApproximateGroundStationPositions( ).at( syntheticDopplerStationName );

    std::vector< double > utcTimes;
    std::vector< Eigen::Vector3d > stationPositions;
    for( double utcOffset : utcOffsets )
    {
        utcTimes.push_back( getSyntheticDopplerStartUtc( ) + utcOffset );
        stationPositions.push_back( stationPosition );
    }
    return timeScaleConverter.getCurrentTimes< double >(
            basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcTimes, stationPositions );
}

std::vector< double > convertTdbTimesToUtcOffsets( const std::vector< double >& tdbTimes )
{
    if( tdbTimes.empty( ) )
    {
        return {};
    }

    earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter;
    Eigen::Vector3d stationPosition = simulation_setup::getCombinedApproximateGroundStationPositions( ).at( syntheticDopplerStationName );
    std::vector< double > utcTimes = timeScaleConverter.getCurrentTimesFromSinglePosition< double >(
            basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, tdbTimes, stationPosition );

    std::vector< double > utcOffsets;
    for( double utcTime : utcTimes )
    {
        utcOffsets.push_back( utcTime - getSyntheticDopplerStartUtc( ) );
    }
    return utcOffsets;
}

std::shared_ptr< SingleObservationSet< double, double > > createSyntheticDopplerObservationSet( const std::vector< double >& utcOffsets,
                                                                                                const double integrationTime )
{
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observations;
    for( double utcOffset : utcOffsets )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > observation( 1 );
        observation( 0 ) = utcOffset;
        observations.push_back( observation );
    }

    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Earth", syntheticDopplerStationName );
    linkEnds[ retransmitter ] = LinkEndId( "SyntheticSpacecraft", "" );
    linkEnds[ receiver ] = LinkEndId( "Earth", syntheticDopplerStationName );

    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
            std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillarySettings->setAncillaryDoubleData( doppler_integration_time, integrationTime );

    return std::make_shared< SingleObservationSet< double, double > >( dsn_n_way_averaged_doppler,
                                                                       linkEnds,
                                                                       observations,
                                                                       convertUtcOffsetsToTdbTimes( utcOffsets ),
                                                                       receiver,
                                                                       std::vector< Eigen::VectorXd >( ),
                                                                       nullptr,
                                                                       ancillarySettings );
}

void checkScalarObservations( const std::shared_ptr< SingleObservationSet< double, double > >& observationSet,
                              const std::vector< double >& expectedObservations )
{
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observations = observationSet->getObservationsReference( );
    BOOST_CHECK_EQUAL( observations.size( ), expectedObservations.size( ) );
    for( unsigned int i = 0; i < expectedObservations.size( ); i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( observations.at( i )( 0 ), expectedObservations.at( i ), 1.0E-14 );
    }
}

void checkUtcOffsets( const std::shared_ptr< SingleObservationSet< double, double > >& observationSet,
                      const std::vector< double >& expectedUtcOffsets )
{
    std::vector< double > utcOffsets = convertTdbTimesToUtcOffsets( observationSet->getObservationTimesReference( ) );
    BOOST_CHECK_EQUAL( utcOffsets.size( ), expectedUtcOffsets.size( ) );
    for( unsigned int i = 0; i < expectedUtcOffsets.size( ); i++ )
    {
        BOOST_CHECK_SMALL( utcOffsets.at( i ) - expectedUtcOffsets.at( i ), 1.0E-6 );
    }
}

}  // namespace

BOOST_AUTO_TEST_CASE( testProcessOdfData )
{
    // Define ODF data paths
    std::string odfFile = tudat::paths::getTudatTestDataPath( ) + "/odf07155.odf";

    // Process ODF file data
    std::string spacecraftName = "MESSENGER";
    auto [ trackingDataList, supplementaryData ] = loadOdfFile( std::vector< std::string >{ odfFile }, spacecraftName, "Earth" );

    spice_interface::loadStandardSpiceKernels( );
    // Create settings for default bodies
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun" };
    std::string globalFrameOrigin = "SSB";
    std::string globalFrameOrientation = "J2000";
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, globalFrameOrigin, globalFrameOrientation );

    // Add high-accuracy Earth settings
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< double, Time >( bodySettings );

    // Create data structure that handles Observed Data in Tudat
    auto observedObservationCollection = observation_models::createObservationCollection< double, Time >( trackingDataList, bodies );

    auto observationSets = observedObservationCollection->getObservationsSets( );

    // Check the observations for NWayRange and
    // DsnNWayAveragedDoppler
    for( const auto& observableTypeEntry : observationSets )
    {
        observation_models::ObservableType observableType = observableTypeEntry.first;
        const auto& linkEndsMap = observableTypeEntry.second;

        for( const auto& linkEndsEntry : linkEndsMap )
        {
            const observation_models::LinkEnds& linkEnds = linkEndsEntry.first;
            const auto& observationSetVector = linkEndsEntry.second;

            for( const auto& observationSet : observationSetVector )
            {
                // Get the observations and times
                auto observations = observationSet->getObservations( );
                auto observationTimes = observationSet->getObservationTimes( );
                auto ancillarySettings = observationSet->getAncillarySettings( );

                if( !observations.empty( ) && !observationTimes.empty( ) )
                {
                    // Check NWayRange observable
                    if( observableType == observation_models::n_way_range )
                    {
                        BOOST_CHECK_EQUAL( observation_models::getObservableName( observableType ), "NWayRange" );
                        BOOST_CHECK_EQUAL( linkEnds.at( transmitter ).getReferencePointName( ), "DSS-14" );
                        BOOST_CHECK_EQUAL( linkEnds.at( retransmitter ).bodyName_, "MESSENGER" );
                        BOOST_CHECK_EQUAL( linkEnds.at( receiver ).getReferencePointName( ), "DSS-14" );
                        BOOST_CHECK_CLOSE_FRACTION( double( observationTimes.at( 0 ) ), 234262616.184812933, 1e-9 );
                        BOOST_CHECK_CLOSE_FRACTION( observations[ 0 ].transpose( )( 0 ), 333589.366953747, 1e-9 );

                        // Check ancillary settings
                        if( ancillarySettings != nullptr )
                        {
                            BOOST_CHECK_EQUAL(
                                    ancillarySettings->getDoubleData( ).at( observation_models::reception_reference_frequency_band ),
                                    1.000000000 );
                            BOOST_CHECK_CLOSE_FRACTION(
                                    ancillarySettings->getDoubleData( ).at( observation_models::sequential_range_lowest_ranging_component ),
                                    14.000000000,
                                    1e-9 );
                            BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::link_ends_delays ).size( ),
                                               3 );
                            BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::frequency_bands ).size( ),
                                               2 );
                        }
                    }

                    // Check DsnNWayAveragedDoppler observable
                    if( observableType == observation_models::dsn_n_way_averaged_doppler )
                    {
                        BOOST_CHECK_EQUAL( observation_models::getObservableName( observableType ), "DsnNWayAveragedDoppler" );
                        BOOST_CHECK_EQUAL( linkEnds.at( transmitter ).getReferencePointName( ), "DSS-14" );
                        BOOST_CHECK_EQUAL( linkEnds.at( retransmitter ).bodyName_, "MESSENGER" );
                        BOOST_CHECK_EQUAL( linkEnds.at( receiver ).getReferencePointName( ), "DSS-14" );
                        BOOST_CHECK_CLOSE_FRACTION( double( observationTimes.at( 0 ) ), 234262457.184812993, 1e-9 );
                        BOOST_CHECK_CLOSE_FRACTION( observations[ 0 ].transpose( )( 0 ), 1.563486099, 1e-9 );

                        // Check ancillary settings
                        if( ancillarySettings != nullptr )
                        {
                            BOOST_CHECK_EQUAL(
                                    ancillarySettings->getDoubleData( ).at( observation_models::reception_reference_frequency_band ),
                                    1.000000000 );
                            BOOST_CHECK_CLOSE_FRACTION(
                                    ancillarySettings->getDoubleData( ).at( observation_models::doppler_integration_time ),
                                    60.000000000,
                                    1e-9 );
                            BOOST_CHECK_CLOSE_FRACTION(
                                    ancillarySettings->getDoubleData( ).at( observation_models::doppler_reference_frequency ),
                                    7177641534.000000000,
                                    1e-9 );
                            BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::link_ends_delays ).size( ),
                                               3 );
                            BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::frequency_bands ).size( ),
                                               2 );
                        }
                    }
                }
                break;
            }
            break;
        }
    }
}

BOOST_AUTO_TEST_CASE( testOdfRampTablesPerFrequencyBand )
{
    // Uplink frequency bands are identified from the DSN uplink frequency allocations
    BOOST_CHECK_EQUAL( getDsnUplinkFrequencyBandName( 2.115E9 ), "S-band" );
    BOOST_CHECK_EQUAL( getDsnUplinkFrequencyBandName( 7.18E9 ), "X-band" );
    BOOST_CHECK_EQUAL( getDsnUplinkFrequencyBandName( 34.4E9 ), "Ka-band" );
    BOOST_CHECK_EQUAL( getDsnUplinkFrequencyBandName( 8.4E9 ), "" );

    // Ramp tables are stored per ground station and uplink band. The file only contains X-band uplinks.
    const std::vector odfFiles{ tudat::paths::getTudatTestDataPath( ) + "/odf07155.odf" };
    std::shared_ptr< input_output::PreProcessedOdfFileContents< Time > > processedOdfFile =
            input_output::processOdfData< Time >( odfFiles, "MESSENGER", false );

    std::set< std::pair< std::string, std::string > > rampTableIds;
    for( auto const& [ rampTableId, rampTable ] : processedOdfFile->getRampTables( ) )
    {
        rampTableIds.insert( rampTableId );
        BOOST_CHECK_EQUAL( rampTable->getFrequencyBand( ), rampTableId.second );
    }
    BOOST_CHECK( rampTableIds ==
                 ( std::set< std::pair< std::string, std::string > >{
                         { "DSS-14", "X-band" }, { "DSS-43", "X-band" }, { "DSS-63", "X-band" } } ) );

    // Supplementary data hold the band-specific ramp tables of each station
    std::vector< std::shared_ptr< data::TrackingSupplementaryData > > supplementaryData =
            loadOdfFile( odfFiles, "MESSENGER", "Earth", false ).second;
    BOOST_CHECK_EQUAL( supplementaryData.size( ), 3u );

    SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    for( const std::shared_ptr< data::TrackingSupplementaryData >& stationSupplementaryData : supplementaryData )
    {
        BOOST_REQUIRE_EQUAL( stationSupplementaryData->getFrequencySupplementaryData( ).size( ), 1u );
        BOOST_CHECK_EQUAL( stationSupplementaryData->getFrequencySupplementaryData( ).at( 0 )->getFrequencyBand( ), "X-band" );
        createGroundStation(
                bodies.at( "Earth" ), stationSupplementaryData->getReferencePointName( ), Eigen::Vector3d( 6378137.0, 0.0, 0.0 ) );
    }

    // Stations must receive a calculator for the X-band only
    setTrackingSupplementaryDataInBodies( bodies, supplementaryData );
    for( const std::string stationName : { "DSS-14", "DSS-43", "DSS-63" } )
    {
        std::shared_ptr< ground_stations::GroundStation > groundStation = bodies.at( "Earth" )->getGroundStation( stationName );
        BOOST_CHECK( groundStation->hasFrequencyCalculator( basic_astrodynamics::x_band ) );
        BOOST_CHECK( !groundStation->hasFrequencyCalculator( basic_astrodynamics::s_band ) );
    }
}

BOOST_AUTO_TEST_CASE( testCompressDopplerDataUsesCadenceRuns )
{
    std::shared_ptr< SingleObservationSet< double, double > > compressedObservationSet = compressDopplerData(
            createSyntheticDopplerObservationSet( { 0.0, 10.0, 20.0, 30.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0 }, 10.0 ), 3 );

    checkScalarObservations( compressedObservationSet, { 10.0, 70.0, 100.0 } );
    checkUtcOffsets( compressedObservationSet, { 10.0, 70.0, 100.0 } );
    BOOST_CHECK_CLOSE_FRACTION(
            compressedObservationSet->getAncillarySettings( )->getAncillaryDoubleData( doppler_integration_time ), 30.0, 1.0E-14 );

    std::shared_ptr< SingleObservationSet< double, double > > gapFreeCompressedObservationSet =
            compressDopplerData( createSyntheticDopplerObservationSet( { 0.0, 10.0, 20.0, 30.0, 40.0, 50.0 }, 10.0 ), 3 );
    checkScalarObservations( gapFreeCompressedObservationSet, { 10.0, 40.0 } );
    checkUtcOffsets( gapFreeCompressedObservationSet, { 10.0, 40.0 } );

    std::shared_ptr< SingleObservationSet< double, double > > shortRunCompressedObservationSet =
            compressDopplerData( createSyntheticDopplerObservationSet( { 0.0, 10.0 }, 10.0 ), 3 );
    BOOST_CHECK_EQUAL( shortRunCompressedObservationSet->getObservationsReference( ).size( ), 0 );
    BOOST_CHECK_EQUAL( shortRunCompressedObservationSet->getObservationTimesReference( ).size( ), 0 );
}

BOOST_AUTO_TEST_SUITE_END( )
