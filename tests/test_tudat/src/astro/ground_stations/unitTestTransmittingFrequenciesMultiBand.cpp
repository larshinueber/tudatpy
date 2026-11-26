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

#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/astro/system_models/vehicleSystems.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ground_stations;

BOOST_AUTO_TEST_SUITE( test_transmitting_frequencies_multi_band )

// DSN uplink frequencies: X-band ramp from an MRO ramp table, and the coherent S-band ramp (221/749)
const double xBandStartFrequency = 7.18322247E9;
const double xBandRampRate = -45.0;
const double sBandStartFrequency = xBandStartFrequency * 221.0 / 749.0;
const double sBandRampRate = xBandRampRate * 221.0 / 749.0;

//! Create a ramp table of two consecutive ramps of 100 s each, starting at the given time.
std::shared_ptr< PiecewiseLinearFrequencyInterpolator > createRampTable( const double startTime,
                                                                         const double startFrequency,
                                                                         const double rampRate )
{
    return std::make_shared< PiecewiseLinearFrequencyInterpolator >(
            std::vector< Time >{ Time( startTime ), Time( startTime + 100.0 ) },
            std::vector< Time >{ Time( startTime + 100.0 ), Time( startTime + 200.0 ) },
            std::vector< double >{ rampRate, 2.0 * rampRate },
            std::vector< double >{ startFrequency, startFrequency + 100.0 * rampRate } );
}

//! Frequency of a ramp table created with createRampTable, at a time within (or after) its second ramp.
double getSecondRampFrequency( const double startTime, const double startFrequency, const double rampRate, const double time )
{
    return startFrequency + 100.0 * rampRate + 2.0 * rampRate * ( time - ( startTime + 100.0 ) );
}

// Second X-band ramp table, to be merged with an X-band table created with createRampTable at t = 0 s
const double secondTableStartTime = 1000.0;
const double secondTableStartFrequency = xBandStartFrequency + 1.0E3;

// Test times within the first (0-200 s) and second (1000-1200 s) ramp table, and the expected frequencies at these times
const double firstTableTestTime = 150.0;
const double secondTableTestTime = 1150.0;
const double xBandFrequencyInFirstTable = getSecondRampFrequency( 0.0, xBandStartFrequency, xBandRampRate, firstTableTestTime );
const double sBandFrequencyInFirstTable = getSecondRampFrequency( 0.0, sBandStartFrequency, sBandRampRate, firstTableTestTime );
// Without merged second table, the frequency at the second test time is extrapolated from the first table
const double unmergedXBandFrequencyInSecondTable = getSecondRampFrequency( 0.0, xBandStartFrequency, xBandRampRate, secondTableTestTime );
const double mergedXBandFrequencyInSecondTable =
        getSecondRampFrequency( secondTableStartTime, secondTableStartFrequency, xBandRampRate, secondTableTestTime );

const double frequencyTolerance = 1.0E-14;

BOOST_AUTO_TEST_CASE( testLegacySingleFrequencyInterpolator )
{
    using namespace tudat::basic_astrodynamics;
    std::shared_ptr< PiecewiseLinearFrequencyInterpolator > xBandInterpolator = createRampTable( 0.0, xBandStartFrequency, xBandRampRate );

    // Band-independent interpolators ignore the frequency band
    ConstantFrequencyInterpolator constantInterpolator( xBandStartFrequency );
    for( const FrequencyBands frequencyBand : { FrequencyBands::s_band, FrequencyBands::x_band, FrequencyBands::ka_band } )
    {
        BOOST_CHECK_EQUAL( ( xBandInterpolator->getTemplatedCurrentFrequency< double, double >( 150.0, frequencyBand ) ),
                           ( xBandInterpolator->getTemplatedCurrentFrequency< double, double >( 150.0 ) ) );
        BOOST_CHECK_EQUAL( ( constantInterpolator.getTemplatedCurrentFrequency< double, double >( 150.0, frequencyBand ) ),
                           xBandStartFrequency );
    }
}

BOOST_AUTO_TEST_CASE( testMultiBandFrequencyInterpolator )
{
    using namespace tudat::basic_astrodynamics;

    std::shared_ptr< PiecewiseLinearFrequencyInterpolator > xBandInterpolator = createRampTable( 0.0, xBandStartFrequency, xBandRampRate );
    std::shared_ptr< PiecewiseLinearFrequencyInterpolator > sBandInterpolator = createRampTable( 0.0, sBandStartFrequency, sBandRampRate );
    MultiBandFrequencyInterpolator multiBandInterpolator( std::map< FrequencyBands, std::shared_ptr< StationFrequencyInterpolator > >{
            { FrequencyBands::x_band, xBandInterpolator }, { FrequencyBands::s_band, sBandInterpolator } } );

    // Frequency and frequency integral of each band are computed by the interpolator of that band
    for( auto const& [ frequencyBand, bandRamps ] : std::map< FrequencyBands, std::shared_ptr< PiecewiseLinearFrequencyInterpolator > >{
                 { FrequencyBands::x_band, xBandInterpolator }, { FrequencyBands::s_band, sBandInterpolator } } )
    {
        for( const double time : { 50.0, 150.0 } )
        {
            BOOST_CHECK_EQUAL( ( multiBandInterpolator.getTemplatedCurrentFrequency< double, double >( time, frequencyBand ) ),
                               ( bandRamps->getTemplatedCurrentFrequency< double, double >( time ) ) );
            BOOST_CHECK_EQUAL( ( multiBandInterpolator.getTemplatedCurrentFrequency< long double, Time >( Time( time ), frequencyBand ) ),
                               ( bandRamps->getTemplatedCurrentFrequency< long double, Time >( Time( time ) ) ) );
        }
        BOOST_CHECK_EQUAL( ( multiBandInterpolator.getTemplatedFrequencyIntegral< double, double >( 50.0, 150.0, frequencyBand ) ),
                           ( bandRamps->getTemplatedFrequencyIntegral< double, double >( 50.0, 150.0 ) ) );
        BOOST_CHECK_EQUAL(
                ( multiBandInterpolator.getTemplatedFrequencyIntegral< long double, Time >( Time( 50.0 ), Time( 150.0 ), frequencyBand ) ),
                ( bandRamps->getTemplatedFrequencyIntegral< long double, Time >( Time( 50.0 ), Time( 150.0 ) ) ) );
    }

    // Bands are different, so the band selection matters
    BOOST_CHECK( ( multiBandInterpolator.getTemplatedCurrentFrequency< double, double >( 50.0, FrequencyBands::x_band ) ) !=
                 ( multiBandInterpolator.getTemplatedCurrentFrequency< double, double >( 50.0, FrequencyBands::s_band ) ) );

    // Without frequency band, the frequency is ambiguous for multiple bands
    BOOST_CHECK_THROW( ( multiBandInterpolator.getTemplatedCurrentFrequency< double, double >( 50.0 ) ), std::runtime_error );
    BOOST_CHECK_THROW( ( multiBandInterpolator.getTemplatedCurrentFrequency< double, double >( 50.0, std::nullopt ) ), std::runtime_error );
    BOOST_CHECK_THROW( ( multiBandInterpolator.getTemplatedFrequencyIntegral< double, double >( 50.0, 150.0 ) ), std::runtime_error );

    // Band without interpolator
    BOOST_CHECK_THROW( ( multiBandInterpolator.getTemplatedCurrentFrequency< double, double >( 50.0, FrequencyBands::ka_band ) ),
                       std::runtime_error );

    // Without frequency band, the frequency is unambiguous for a single band
    MultiBandFrequencyInterpolator singleBandMultiFrequencyInterpolator(
            std::map< FrequencyBands, std::shared_ptr< StationFrequencyInterpolator > >{ { FrequencyBands::x_band, xBandInterpolator } } );
    BOOST_CHECK_EQUAL( ( singleBandMultiFrequencyInterpolator.getTemplatedCurrentFrequency< double, double >( 50.0 ) ),
                       ( xBandInterpolator->getTemplatedCurrentFrequency< double, double >( 50.0 ) ) );
    BOOST_CHECK_EQUAL(
            ( singleBandMultiFrequencyInterpolator.getTemplatedFrequencyIntegral< double, double >( 50.0, 150.0, std::nullopt ) ),
            ( xBandInterpolator->getTemplatedFrequencyIntegral< double, double >( 50.0, 150.0 ) ) );
}

BOOST_AUTO_TEST_CASE( testCombineTransmittingFrequencyCalculatorsWithoutExistingCalculator )
{
    using namespace tudat::basic_astrodynamics;

    std::shared_ptr< PiecewiseLinearFrequencyInterpolator > xBandRampedInterpolator =
            createRampTable( 0.0, xBandStartFrequency, xBandRampRate );
    auto sBandConstantInterpolator = std::make_shared< ConstantFrequencyInterpolator >( sBandStartFrequency );
    auto multiBandCalculator =
            std::make_shared< MultiBandFrequencyInterpolator >( std::map< FrequencyBands, std::shared_ptr< StationFrequencyInterpolator > >{
                    { FrequencyBands::x_band, xBandRampedInterpolator } } );

    // Without band: calculator to add is used as is
    {
        std::shared_ptr< StationFrequencyInterpolator > combinedCalculator =
                combineTransmittingFrequencyCalculators( nullptr, xBandRampedInterpolator, std::nullopt );
        BOOST_REQUIRE( combinedCalculator != nullptr );
        BOOST_CHECK( std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >( combinedCalculator ) == nullptr );
        BOOST_CHECK_CLOSE_FRACTION( ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime ) ),
                                    xBandFrequencyInFirstTable,
                                    frequencyTolerance );
    }
    {
        std::shared_ptr< StationFrequencyInterpolator > combinedCalculator =
                combineTransmittingFrequencyCalculators( nullptr, sBandConstantInterpolator, std::nullopt );
        BOOST_REQUIRE( combinedCalculator != nullptr );
        BOOST_CHECK( std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >( combinedCalculator ) == nullptr );
        BOOST_CHECK_EQUAL( ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime ) ),
                           sBandStartFrequency );
    }
    {
        std::shared_ptr< MultiBandFrequencyInterpolator > combinedCalculator = std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >(
                combineTransmittingFrequencyCalculators( nullptr, multiBandCalculator, std::nullopt ) );
        BOOST_REQUIRE( combinedCalculator != nullptr );
        BOOST_CHECK( combinedCalculator->getAvailableBands( ) == std::vector< FrequencyBands >{ FrequencyBands::x_band } );
        BOOST_CHECK_CLOSE_FRACTION(
                ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime, FrequencyBands::x_band ) ),
                xBandFrequencyInFirstTable,
                frequencyTolerance );
    }

    // With band: calculator to add is stored for that band in a new multi-band calculator
    {
        std::shared_ptr< MultiBandFrequencyInterpolator > combinedCalculator = std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >(
                combineTransmittingFrequencyCalculators( nullptr, xBandRampedInterpolator, FrequencyBands::x_band ) );
        BOOST_REQUIRE( combinedCalculator != nullptr );
        BOOST_CHECK( combinedCalculator->getAvailableBands( ) == std::vector< FrequencyBands >{ FrequencyBands::x_band } );
        BOOST_CHECK_CLOSE_FRACTION(
                ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime, FrequencyBands::x_band ) ),
                xBandFrequencyInFirstTable,
                frequencyTolerance );
    }
    {
        std::shared_ptr< MultiBandFrequencyInterpolator > combinedCalculator = std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >(
                combineTransmittingFrequencyCalculators( nullptr, sBandConstantInterpolator, FrequencyBands::s_band ) );
        BOOST_REQUIRE( combinedCalculator != nullptr );
        BOOST_CHECK( combinedCalculator->getAvailableBands( ) == std::vector< FrequencyBands >{ FrequencyBands::s_band } );
        BOOST_CHECK_EQUAL(
                ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime, FrequencyBands::s_band ) ),
                sBandStartFrequency );
    }

    // A multi-band calculator cannot be stored for a single band
    BOOST_CHECK_THROW( combineTransmittingFrequencyCalculators( nullptr, multiBandCalculator, FrequencyBands::x_band ),
                       std::runtime_error );
}

BOOST_AUTO_TEST_CASE( testCombineTransmittingFrequencyCalculatorsWithBandIndependentCalculator )
{
    using namespace tudat::basic_astrodynamics;

    std::shared_ptr< PiecewiseLinearFrequencyInterpolator > secondRampedInterpolator =
            createRampTable( secondTableStartTime, secondTableStartFrequency, xBandRampRate );
    auto constantCalculator = std::make_shared< ConstantFrequencyInterpolator >( xBandStartFrequency );
    auto multiBandCalculator =
            std::make_shared< MultiBandFrequencyInterpolator >( std::map< FrequencyBands, std::shared_ptr< StationFrequencyInterpolator > >{
                    { FrequencyBands::x_band, secondRampedInterpolator } } );

    // Piecewise linear + piecewise linear: ramp tables are merged into the existing calculator, in place
    {
        std::shared_ptr< PiecewiseLinearFrequencyInterpolator > existingRampedFrequencyInterpolator =
                createRampTable( 0.0, xBandStartFrequency, xBandRampRate );
        std::shared_ptr< StationFrequencyInterpolator > previouslyRetrievedCalculator = existingRampedFrequencyInterpolator;

        std::shared_ptr< StationFrequencyInterpolator > combinedCalculator =
                combineTransmittingFrequencyCalculators( existingRampedFrequencyInterpolator, secondRampedInterpolator, std::nullopt );
        BOOST_REQUIRE( combinedCalculator != nullptr );
        BOOST_CHECK( std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >( combinedCalculator ) == nullptr );
        BOOST_CHECK_CLOSE_FRACTION( ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime ) ),
                                    xBandFrequencyInFirstTable,
                                    frequencyTolerance );
        BOOST_CHECK_CLOSE_FRACTION( ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( secondTableTestTime ) ),
                                    mergedXBandFrequencyInSecondTable,
                                    frequencyTolerance );

        // Objects referring to the existing calculator see the merged ramps
        BOOST_CHECK_CLOSE_FRACTION(
                ( previouslyRetrievedCalculator->getTemplatedCurrentFrequency< double, double >( secondTableTestTime ) ),
                mergedXBandFrequencyInSecondTable,
                frequencyTolerance );
    }

    // Combinations that cannot be merged: an exception is thrown, and the existing calculator is unchanged
    {
        std::shared_ptr< PiecewiseLinearFrequencyInterpolator > existingRampedFrequencyInterpolator =
                createRampTable( 0.0, xBandStartFrequency, xBandRampRate );

        // Only piecewise linear calculators can be merged
        BOOST_CHECK_THROW( combineTransmittingFrequencyCalculators( existingRampedFrequencyInterpolator, constantCalculator, std::nullopt ),
                           std::runtime_error );

        // Band-specific calculators cannot be added to a band-independent calculator
        BOOST_CHECK_THROW( combineTransmittingFrequencyCalculators(
                                   existingRampedFrequencyInterpolator, secondRampedInterpolator, FrequencyBands::x_band ),
                           std::runtime_error );
        BOOST_CHECK_THROW(
                combineTransmittingFrequencyCalculators( existingRampedFrequencyInterpolator, constantCalculator, FrequencyBands::x_band ),
                std::runtime_error );
        BOOST_CHECK_THROW(
                combineTransmittingFrequencyCalculators( existingRampedFrequencyInterpolator, multiBandCalculator, std::nullopt ),
                std::runtime_error );

        BOOST_CHECK_CLOSE_FRACTION(
                ( existingRampedFrequencyInterpolator->getTemplatedCurrentFrequency< double, double >( secondTableTestTime ) ),
                unmergedXBandFrequencyInSecondTable,
                frequencyTolerance );
    }

    {
        auto existingConstantCalculator = std::make_shared< ConstantFrequencyInterpolator >( xBandStartFrequency );

        // Only piecewise linear calculators can be merged
        BOOST_CHECK_THROW( combineTransmittingFrequencyCalculators( existingConstantCalculator, secondRampedInterpolator, std::nullopt ),
                           std::runtime_error );
        BOOST_CHECK_THROW( combineTransmittingFrequencyCalculators( existingConstantCalculator, constantCalculator, std::nullopt ),
                           std::runtime_error );

        // Band-specific calculators cannot be added to a band-independent calculator
        BOOST_CHECK_THROW(
                combineTransmittingFrequencyCalculators( existingConstantCalculator, secondRampedInterpolator, FrequencyBands::x_band ),
                std::runtime_error );
        BOOST_CHECK_THROW( combineTransmittingFrequencyCalculators( existingConstantCalculator, multiBandCalculator, std::nullopt ),
                           std::runtime_error );

        BOOST_CHECK_EQUAL( ( existingConstantCalculator->getTemplatedCurrentFrequency< double, double >( secondTableTestTime ) ),
                           xBandStartFrequency );
    }
}

BOOST_AUTO_TEST_CASE( testCombineTransmittingFrequencyCalculatorsWithMultiBandCalculator )
{
    using namespace tudat::basic_astrodynamics;

    // Existing multi-band calculator, with a piecewise linear X-band ramp table
    auto createXbandMultiCalculator = [ & ]( ) {
        return std::make_shared< MultiBandFrequencyInterpolator >(
                std::map< FrequencyBands, std::shared_ptr< StationFrequencyInterpolator > >{
                        { FrequencyBands::x_band, createRampTable( 0.0, xBandStartFrequency, xBandRampRate ) } } );
    };

    // Calculator for a new band is added, and the existing band is unchanged
    {
        std::shared_ptr< MultiBandFrequencyInterpolator > combinedFrequencyInterpolator =
                std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >(
                        combineTransmittingFrequencyCalculators( createXbandMultiCalculator( ),
                                                                 createRampTable( 0.0, sBandStartFrequency, sBandRampRate ),
                                                                 FrequencyBands::s_band ) );
        BOOST_REQUIRE( combinedFrequencyInterpolator != nullptr );
        BOOST_CHECK( combinedFrequencyInterpolator->getAvailableBands( ) ==
                     ( std::vector< FrequencyBands >{ FrequencyBands::s_band, FrequencyBands::x_band } ) );
        BOOST_CHECK_CLOSE_FRACTION( ( combinedFrequencyInterpolator->getTemplatedCurrentFrequency< double, double >(
                                            firstTableTestTime, FrequencyBands::s_band ) ),
                                    sBandFrequencyInFirstTable,
                                    frequencyTolerance );
        BOOST_CHECK_CLOSE_FRACTION( ( combinedFrequencyInterpolator->getTemplatedCurrentFrequency< double, double >(
                                            secondTableTestTime, FrequencyBands::x_band ) ),
                                    unmergedXBandFrequencyInSecondTable,
                                    frequencyTolerance );
    }
    {
        std::shared_ptr< MultiBandFrequencyInterpolator > combinedFrequencyInterpolator =
                std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >(
                        combineTransmittingFrequencyCalculators( createXbandMultiCalculator( ),
                                                                 std::make_shared< ConstantFrequencyInterpolator >( sBandStartFrequency ),
                                                                 FrequencyBands::s_band ) );
        BOOST_REQUIRE( combinedFrequencyInterpolator != nullptr );
        BOOST_CHECK( combinedFrequencyInterpolator->getAvailableBands( ) ==
                     ( std::vector< FrequencyBands >{ FrequencyBands::s_band, FrequencyBands::x_band } ) );
        BOOST_CHECK_EQUAL( ( combinedFrequencyInterpolator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime,
                                                                                                            FrequencyBands::s_band ) ),
                           sBandStartFrequency );
        BOOST_CHECK_CLOSE_FRACTION( ( combinedFrequencyInterpolator->getTemplatedCurrentFrequency< double, double >(
                                            secondTableTestTime, FrequencyBands::x_band ) ),
                                    unmergedXBandFrequencyInSecondTable,
                                    frequencyTolerance );
    }

    // Piecewise linear calculator for an existing piecewise linear band is merged into it, in place
    {
        std::shared_ptr< MultiBandFrequencyInterpolator > existingCalculator = createXbandMultiCalculator( );
        std::shared_ptr< StationFrequencyInterpolator > previouslyRetrievedCalculator =
                existingCalculator->getBandInterpolator( FrequencyBands::x_band );

        std::shared_ptr< MultiBandFrequencyInterpolator > combinedCalculator = std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >(
                combineTransmittingFrequencyCalculators( existingCalculator,
                                                         createRampTable( secondTableStartTime, secondTableStartFrequency, xBandRampRate ),
                                                         FrequencyBands::x_band ) );
        BOOST_REQUIRE( combinedCalculator != nullptr );
        BOOST_CHECK( combinedCalculator->getAvailableBands( ) == std::vector< FrequencyBands >{ FrequencyBands::x_band } );
        BOOST_CHECK_CLOSE_FRACTION(
                ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime, FrequencyBands::x_band ) ),
                xBandFrequencyInFirstTable,
                frequencyTolerance );
        BOOST_CHECK_CLOSE_FRACTION(
                ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( secondTableTestTime, FrequencyBands::x_band ) ),
                mergedXBandFrequencyInSecondTable,
                frequencyTolerance );

        // Objects referring to the existing calculator see the merged ramps
        BOOST_CHECK_CLOSE_FRACTION(
                ( previouslyRetrievedCalculator->getTemplatedCurrentFrequency< double, double >( secondTableTestTime ) ),
                mergedXBandFrequencyInSecondTable,
                frequencyTolerance );
    }

    // Multi-band calculator is merged band by band: new bands are added, existing piecewise linear bands are merged
    {
        auto calculatorToAdd = std::make_shared< MultiBandFrequencyInterpolator >(
                std::map< FrequencyBands, std::shared_ptr< StationFrequencyInterpolator > >{
                        { FrequencyBands::s_band, std::make_shared< ConstantFrequencyInterpolator >( sBandStartFrequency ) },
                        { FrequencyBands::x_band, createRampTable( secondTableStartTime, secondTableStartFrequency, xBandRampRate ) } } );

        std::shared_ptr< MultiBandFrequencyInterpolator > combinedCalculator = std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >(
                combineTransmittingFrequencyCalculators( createXbandMultiCalculator( ), calculatorToAdd, std::nullopt ) );
        BOOST_REQUIRE( combinedCalculator != nullptr );
        BOOST_CHECK( combinedCalculator->getAvailableBands( ) ==
                     ( std::vector< FrequencyBands >{ FrequencyBands::s_band, FrequencyBands::x_band } ) );
        BOOST_CHECK_EQUAL(
                ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime, FrequencyBands::s_band ) ),
                sBandStartFrequency );
        BOOST_CHECK_CLOSE_FRACTION(
                ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime, FrequencyBands::x_band ) ),
                xBandFrequencyInFirstTable,
                frequencyTolerance );
        BOOST_CHECK_CLOSE_FRACTION(
                ( combinedCalculator->getTemplatedCurrentFrequency< double, double >( secondTableTestTime, FrequencyBands::x_band ) ),
                mergedXBandFrequencyInSecondTable,
                frequencyTolerance );
    }

    // Combinations that cannot be merged: an exception is thrown, and the existing calculator is unchanged
    {
        std::shared_ptr< MultiBandFrequencyInterpolator > existingCalculator = createXbandMultiCalculator( );

        // Band-independent calculators cannot be added to a multi-band calculator
        BOOST_CHECK_THROW(
                combineTransmittingFrequencyCalculators( existingCalculator,
                                                         createRampTable( secondTableStartTime, secondTableStartFrequency, xBandRampRate ),
                                                         std::nullopt ),
                std::runtime_error );
        BOOST_CHECK_THROW(
                combineTransmittingFrequencyCalculators(
                        existingCalculator, std::make_shared< ConstantFrequencyInterpolator >( xBandStartFrequency ), std::nullopt ),
                std::runtime_error );

        // Only piecewise linear calculators can be merged into an existing band
        BOOST_CHECK_THROW(
                combineTransmittingFrequencyCalculators( existingCalculator,
                                                         std::make_shared< ConstantFrequencyInterpolator >( xBandStartFrequency ),
                                                         FrequencyBands::x_band ),
                std::runtime_error );

        // A multi-band calculator cannot be stored for a single band
        auto multiBandCalculatorToAdd = std::make_shared< MultiBandFrequencyInterpolator >(
                std::map< FrequencyBands, std::shared_ptr< StationFrequencyInterpolator > >{
                        { FrequencyBands::s_band, std::make_shared< ConstantFrequencyInterpolator >( sBandStartFrequency ) } } );
        BOOST_CHECK_THROW( combineTransmittingFrequencyCalculators( existingCalculator, multiBandCalculatorToAdd, FrequencyBands::s_band ),
                           std::runtime_error );

        // Multi-band calculator with a band that cannot be merged: none of its bands is added
        auto partiallyInvalidCalculatorToAdd = std::make_shared< MultiBandFrequencyInterpolator >(
                std::map< FrequencyBands, std::shared_ptr< StationFrequencyInterpolator > >{
                        { FrequencyBands::s_band, std::make_shared< ConstantFrequencyInterpolator >( sBandStartFrequency ) },
                        { FrequencyBands::x_band, std::make_shared< ConstantFrequencyInterpolator >( xBandStartFrequency ) } } );
        BOOST_CHECK_THROW( combineTransmittingFrequencyCalculators( existingCalculator, partiallyInvalidCalculatorToAdd, std::nullopt ),
                           std::runtime_error );

        // Null calculator cannot be added
        BOOST_CHECK_THROW( combineTransmittingFrequencyCalculators( existingCalculator, nullptr, FrequencyBands::s_band ),
                           std::runtime_error );

        BOOST_CHECK( existingCalculator->getAvailableBands( ) == std::vector< FrequencyBands >{ FrequencyBands::x_band } );
        BOOST_CHECK_CLOSE_FRACTION(
                ( existingCalculator->getTemplatedCurrentFrequency< double, double >( secondTableTestTime, FrequencyBands::x_band ) ),
                unmergedXBandFrequencyInSecondTable,
                frequencyTolerance );
    }

    // Piecewise linear calculator cannot be merged into an existing band with a constant calculator
    {
        auto existingCalculator = std::make_shared< MultiBandFrequencyInterpolator >(
                std::map< FrequencyBands, std::shared_ptr< StationFrequencyInterpolator > >{
                        { FrequencyBands::x_band, std::make_shared< ConstantFrequencyInterpolator >( xBandStartFrequency ) } } );

        BOOST_CHECK_THROW(
                combineTransmittingFrequencyCalculators( existingCalculator,
                                                         createRampTable( secondTableStartTime, secondTableStartFrequency, xBandRampRate ),
                                                         FrequencyBands::x_band ),
                std::runtime_error );
        BOOST_CHECK( existingCalculator->getAvailableBands( ) == std::vector< FrequencyBands >{ FrequencyBands::x_band } );
        BOOST_CHECK_EQUAL(
                ( existingCalculator->getTemplatedCurrentFrequency< double, double >( secondTableTestTime, FrequencyBands::x_band ) ),
                xBandStartFrequency );
    }
}

BOOST_AUTO_TEST_CASE( testGroundStationAndVehicleFrequencyCalculatorAccess )
{
    using namespace tudat::basic_astrodynamics;

    std::shared_ptr< PiecewiseLinearFrequencyInterpolator > xBandRampedFrequencyInterpolator =
            createRampTable( 0.0, xBandStartFrequency, xBandRampRate );
    auto sBandCalculator = std::make_shared< ConstantFrequencyInterpolator >( sBandStartFrequency );

    // Ground station with band-specific calculators
    GroundStation multiBandStation( std::make_shared< GroundStationState >( Eigen::Vector3d( 6378137.0, 0.0, 0.0 ) ), nullptr, "DSS-43" );
    BOOST_CHECK( !multiBandStation.hasFrequencyCalculator( ) );
    BOOST_CHECK( !multiBandStation.hasFrequencyCalculator( FrequencyBands::x_band ) );

    multiBandStation.addTransmittingFrequencyCalculator( xBandRampedFrequencyInterpolator, FrequencyBands::x_band );
    multiBandStation.addTransmittingFrequencyCalculator( sBandCalculator, FrequencyBands::s_band );
    BOOST_CHECK( std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >( multiBandStation.getTransmittingFrequencyCalculator( ) ) !=
                 nullptr );

    // Ground station has S/X-band calculators, but not Ka-band
    BOOST_CHECK( multiBandStation.hasFrequencyCalculator( ) );
    BOOST_CHECK( multiBandStation.hasFrequencyCalculator( FrequencyBands::x_band ) );
    BOOST_CHECK( multiBandStation.hasFrequencyCalculator( FrequencyBands::s_band ) );
    BOOST_CHECK( !multiBandStation.hasFrequencyCalculator( FrequencyBands::ka_band ) );

    // Frequency of each band is computed by the calculator of that band
    BOOST_CHECK_CLOSE_FRACTION( ( multiBandStation.getTransmittingFrequencyCalculator( FrequencyBands::x_band )
                                          ->getTemplatedCurrentFrequency< double, double >( firstTableTestTime ) ),
                                xBandFrequencyInFirstTable,
                                frequencyTolerance );
    BOOST_CHECK_EQUAL( ( multiBandStation.getTransmittingFrequencyCalculator( FrequencyBands::s_band )
                                 ->getTemplatedCurrentFrequency< double, double >( firstTableTestTime ) ),
                       sBandStartFrequency );

    // Retrieving a calculator for a band without a calculator, or adding a calculator for an existing band, throws an exception
    BOOST_CHECK_THROW( multiBandStation.getTransmittingFrequencyCalculator( FrequencyBands::ka_band ), std::runtime_error );
    BOOST_CHECK_THROW( multiBandStation.addTransmittingFrequencyCalculator( sBandCalculator ), std::runtime_error );

    // Ground station with a band-independent calculator: calculator is used for any band
    GroundStation singleBandStation( std::make_shared< GroundStationState >( Eigen::Vector3d( 6378137.0, 0.0, 0.0 ) ), nullptr, "DSS-14" );
    singleBandStation.addTransmittingFrequencyCalculator( xBandRampedFrequencyInterpolator );
    BOOST_CHECK( std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >( singleBandStation.getTransmittingFrequencyCalculator( ) ) ==
                 nullptr );
    for( const FrequencyBands frequencyBand : { FrequencyBands::s_band, FrequencyBands::x_band, FrequencyBands::ka_band } )
    {
        BOOST_CHECK( singleBandStation.hasFrequencyCalculator( frequencyBand ) );
        BOOST_CHECK_CLOSE_FRACTION( ( singleBandStation.getTransmittingFrequencyCalculator( frequencyBand )
                                              ->getTemplatedCurrentFrequency< double, double >( firstTableTestTime ) ),
                                    xBandFrequencyInFirstTable,
                                    frequencyTolerance );
    }

    // Vehicle systems with band-specific calculators
    system_models::VehicleSystems vehicleSystems;
    vehicleSystems.addTransmittedFrequencyCalculator( xBandRampedFrequencyInterpolator, FrequencyBands::x_band );
    BOOST_CHECK_THROW( vehicleSystems.addTransmittedFrequencyCalculator( sBandCalculator ), std::runtime_error );
    std::shared_ptr< MultiBandFrequencyInterpolator > vehicleCalculator =
            std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >( vehicleSystems.getTransmittedFrequencyCalculator( ) );
    BOOST_REQUIRE( vehicleCalculator != nullptr );
    BOOST_CHECK( vehicleCalculator->getAvailableBands( ) == std::vector< FrequencyBands >{ FrequencyBands::x_band } );
    BOOST_CHECK_CLOSE_FRACTION(
            ( vehicleCalculator->getTemplatedCurrentFrequency< double, double >( firstTableTestTime, FrequencyBands::x_band ) ),
            xBandFrequencyInFirstTable,
            frequencyTolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
