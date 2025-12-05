/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "GIDI.hpp"

namespace GIDI {

namespace Transporting {

/*! \class Settings
 * This class is used to instruct deterministic methods and the Monte Carlo API MCGIDI on what data are being requested.
*/

/* *********************************************************************************************************//**
 * @param a_projectileID            [in]    The PoPs id for the projectile.
 * @param a_delayedNeutrons         [in]    Flag indicating whether or not delayed neutron data are to be include in the requested data.
 ***********************************************************************************************************/

Settings::Settings( std::string const &a_projectileID, DelayedNeutrons a_delayedNeutrons ) :
        m_projectileID( a_projectileID ),
        m_delayedNeutrons( a_delayedNeutrons ),
        m_nuclearPlusCoulombInterferenceOnly( false ),
        m_throwOnError( true ),
        m_zeroDepositionIfAllProductsTracked( true ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Settings::~Settings( ) {

}

/* *********************************************************************************************************//**
 * Returns a Vector of 0.0's of the proper length for the projectile's multi-group data.
 *
 * @param a_particles   [in]    The list of particles to be transported.
 * @param a_collapse    [in]    If true, the length of the returned vector is for the collapsed multi-group, otherwise, it is for the uncollapsed multi-group.
 *
 * @return                      The Vector of 0.0.
 ***********************************************************************************************************/

Vector Settings::multiGroupZeroVector( Particles const &a_particles, bool a_collapse ) const {

    Particle const *projectile( a_particles.particle( projectileID( ) ) );

    std::size_t n1 = projectile->fineMultiGroup( ).numberOfGroups( );
    if( a_collapse ) n1 = projectile->numberOfGroups( );

    Vector vector( n1 );
    return( vector );
}

/* *********************************************************************************************************//**
 * Returns a Matrix of 0.0's of the proper length for the projectile's and product's multi-group data.
 *
 * @param a_particles               [in]    The list of particles to be transported.
 * @param a_pid                     [in]    The PoPs index for the product.
 * @param a_collapse                [in]    If true, the length of the returned vector is for the collapsed multi-group, otherwise, it is for the uncollapsed multi-group.
 *
 * @return                                  The Matrix of 0.0.
 ***********************************************************************************************************/

Matrix Settings::multiGroupZeroMatrix( Particles const &a_particles, std::string const &a_pid, bool a_collapse ) const {

    Particle const *projectile( a_particles.particle( projectileID( ) ) );
    Particle const *product( a_particles.particle( a_pid ) );

    std::size_t n1 = projectile->fineMultiGroup( ).numberOfGroups( );
    std::size_t n2 = product->fineMultiGroup( ).numberOfGroups( );
    if( a_collapse ) {
        n1 = projectile->numberOfGroups( );
        n2 = product->numberOfGroups( );
    }

    Matrix matrix( n1, n2 );
    return( matrix );
}

#if 0
/* *********************************************************************************************************//**
 * Prints the contents of *this* to std::cout. Mainly used for debugging.
 ***********************************************************************************************************/

void Settings::print( ) const {

    std::cout << "setting info:" << std::endl;
    std::cout << "  transport mode = " << m_mode << std::endl;

    std::cout << "  delayed neutrons ";
    if( m_delayedNeutrons ) {
        std::cout << "on" << std::endl; }
    else {
        std::cout << "off" << std::endl;
    }

    std::cout << "  throw on error ";
    if( m_throwOnError ) {
        std::cout << "on" << std::endl; }
    else {
        std::cout << "off" << std::endl;
    }
}
#endif

/*! \class MG
 * This class is used to instruct deterministic methods on what data are being requested.
*/

/* *********************************************************************************************************//**
 * @param a_projectileID            [in]    The PoPs index for the projectile.
 * @param a_mode                    [in]    Specifies the type of data to use or retrieve for transport codes.
 * @param a_delayedNeutrons         [in]    Flag indicating whether or not delayed neutron data are to be include in the requested data.
 ***********************************************************************************************************/

MG::MG( std::string const &a_projectileID, Mode a_mode, DelayedNeutrons a_delayedNeutrons ) :
        Settings( a_projectileID, a_delayedNeutrons ),
        m_mode( a_mode ),
        m_useMultiGroupSummedData( true ) {

}

/* *********************************************************************************************************//**
 * Searches the suite *a_suite* for the form style specified by *m_mode* and matching one in *a_temperatureInfo*.
 * This only works for multi-group data (i.e., multiGroup or multiGroupWithSnElasticUpScatter type data).
 * If *a_label* is not an empty string, then it is used as the label for the form to return and *a_temperatureInfo* is ignored.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_suite               [in]    The suite to search for the requested form.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_dataType            [in]    The type of data being required. Only used if data not found.
 * @param a_label               [in]    The not an empty string, this is used as the label for the form to return.
 ***********************************************************************************************************/

Form const *MG::form( LUPI::StatusMessageReporting &a_smr, GIDI::Suite const &a_suite, Styles::TemperatureInfo const &a_temperatureInfo,
                    std::string a_dataType, std::string const &a_label ) const {

    std::string label;

    if( a_label != "" ) {
        label = a_label; }
    else if( m_mode == Mode::multiGroup ) {
        label = a_temperatureInfo.heatedMultiGroup( ); }
    else if( m_mode == Mode::multiGroupWithSnElasticUpScatter ) {
        label = a_temperatureInfo.SnElasticUpScatter( );
    }

    Suite::const_iterator iter = a_suite.find( label, true );
    if( ( iter == a_suite.end( ) ) && ( a_label == "" ) ) {
        if( m_mode == Mode::multiGroupWithSnElasticUpScatter ) iter = a_suite.find( a_temperatureInfo.heatedMultiGroup( ), true );
    }

    if( iter == a_suite.end( ) ) {
        std::string sourceInfo( a_suite.toXLink( ) );
        ProtareSingle const *rootProtareSingle = static_cast<ProtareSingle const *>( a_suite.root( ) );
        sourceInfo += " for protare " + rootProtareSingle->projectile( ).ID( ) + " + " + rootProtareSingle->target( ).ID( ) 
                + " for " + rootProtareSingle->evaluation( ) + ".";

        if( throwOnError( ) ) {
            throw Exception( "ERROR from GIDI::MG::form: label '" + label + "' not found in suite '" + sourceInfo ); }
        else {
            std::string warning( "data for " + a_dataType + " not found with label '" + label + "' in suite " + sourceInfo );
            smr_setReportError2p( a_smr.smr( ), 0, 0, warning.c_str( ) );
        }
        return( nullptr );
    }

    return( *iter );
}

}

}
