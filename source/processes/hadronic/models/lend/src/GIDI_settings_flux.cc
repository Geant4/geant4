/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "GIDI.hpp"

#include <ptwX.h>
#include <ptwXY.h>

namespace GIDI {

namespace Transporting {

/*! \class Flux_order
 * Specifies the flux data for a specified Legendre order (see class Flux).
 */

/* *********************************************************************************************************//**
 * @param a_order           [in]    The Legendre order of the flux data.
 * @param a_length          [in]    The number of *a_energies* values.
 * @param a_energies        [in]    The list of energies that the flux *a_fluxes* are given at.
 * @param a_fluxes          [in]    The flux for this Legendre order.
 ***********************************************************************************************************/

Flux_order::Flux_order( int a_order, int a_length, double const *a_energies, double const *a_fluxes ) :
        m_order( a_order ) {

    for( int i1 = 0; i1 < a_length; ++i1 ) m_energies.push_back( a_energies[i1] );
    for( int i1 = 0; i1 < a_length; ++i1 ) m_fluxes.push_back( a_fluxes[i1] );
}

/* *********************************************************************************************************//**
 * @param a_order           [in]    The Legendre order of the flux data.
 * @param a_energies        [in]    The list of energies that the flux *a_fluxes* are given at.
 * @param a_fluxes          [in]    The flux for this Legendre order.
 ***********************************************************************************************************/

Flux_order::Flux_order( int a_order, std::vector<double> const &a_energies, std::vector<double> const &a_fluxes ) :
        m_order( a_order ),
        m_energies( a_energies ),
        m_fluxes( a_fluxes ) {

    if( a_energies.size( ) != a_fluxes.size( ) ) throw Exception( "Flux_order::Flux_order: a_energies.size( ) != a_fluxes.size( )." );
}

/* *********************************************************************************************************//**
 * @param a_fluxOrder       [in]    The flux order to copy.
 ***********************************************************************************************************/

Flux_order::Flux_order( Flux_order const &a_fluxOrder ) :
        m_order( a_fluxOrder.order( ) ),
        m_energies( a_fluxOrder.v_energies( ) ),
        m_fluxes( a_fluxOrder.v_fluxes( ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Flux_order::~Flux_order( ) {

}

/* *********************************************************************************************************//**
 * Print the Flux_order to std::cout. Mainly for debugging.
 *
 * @param a_valuesPerLine         [in]    The number of points (i.e., energy, flux pairs) to print per line.
 ***********************************************************************************************************/

void Flux_order::print( int a_valuesPerLine ) const {

    int nE = (int) m_energies.size( );
    bool printIndent = true;

    std::cout << "    ORDER: " << m_order << "  (number of points = " << m_energies.size( ) << ")" << std::endl;
    for( int iE = 0; iE < nE; ++iE ) {
        if( printIndent ) std::cout << "    ";
        printIndent = false;
        std::string buffer = LUPI::Misc::argumentsToString( "   %15.8e %15.8e", m_energies[iE], m_fluxes[iE] );
        std::cout << buffer;
        if( ( ( iE + 1 ) % a_valuesPerLine ) == 0 ) {
            std::cout << std::endl;
            printIndent = true;
        }
    }
    if( nE % a_valuesPerLine ) std::cout << std::endl;
}

/*! \class Flux
 * Specifies the flux data as a list of Flux_order's.
 */

/* *********************************************************************************************************//**
 * @param a_label               [in]    The label for the flux.
 * @param a_temperature         [in]    The temperature the 
 ***********************************************************************************************************/

Flux::Flux( std::string const &a_label, double a_temperature ) :
        m_label( a_label ),
        m_temperature( a_temperature ) {

}

/* *********************************************************************************************************//**
 * @param a_label               [in]    The label for the flux.
 * @param a_temperature         [in]    The temperature the 
 ***********************************************************************************************************/

Flux::Flux( char const *a_label, double a_temperature ) :
        m_label( a_label ),
        m_temperature( a_temperature ) {

}

/* *********************************************************************************************************//**
 * @param a_flux                [in]    The Flux to copy.
 ***********************************************************************************************************/

Flux::Flux( Flux const &a_flux ) :
        m_label( a_flux.label( ) ),
        m_temperature( a_flux.temperature( ) ) {

    for( int i1 = 0; i1 <= a_flux.maxOrder( ); ++i1 ) { addFluxOrder( a_flux[i1] ); }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Flux::~Flux( ) {

}

/* *********************************************************************************************************//**
 * Adds a Flux_order. The Flux_order's must be added sequentially.
 *
 * @param a_fluxOrder           [in]    The Flux_order to add to *this*.
 ***********************************************************************************************************/

void Flux::addFluxOrder( Flux_order const &a_fluxOrder ) {
/*
*   Orders can only be added in sequence (e.g., 0 first, then 1, ...).
*/
    int order = a_fluxOrder.order( );

    if( order > (int) m_fluxOrders.size( ) ) throw Exception( "Flux::addFluxOrder: order > (int) m_fluxOrders.size( )." );
    m_fluxOrders.push_back( a_fluxOrder );
}

/* *********************************************************************************************************//**
 * Multi-groups the flux and returns the result.
 *
 * @param a_multiGroup          [in]    The Flux to copy.
 * @return                      [in]    The Multi-group flux.
 ***********************************************************************************************************/

ProcessedFlux Flux::process( std::vector<double> const &a_multiGroup ) const {
/*
    Currently only does l=0 flux.
*/
    int i1 = 0;
    std::vector<double> groupedFlux;

    ptwXPoints *boundaries = ptwX_create( nullptr, a_multiGroup.size( ), a_multiGroup.size( ), &(a_multiGroup[0]) );
    if( boundaries == nullptr ) throw Exception( "ptwX_create failted for boundaries." );

    for( ; i1 < 1; ++i1 ) {     // only do l=0 currenlty hence ' i1 < 1' test.
        Flux_order const *__fluxOrder = &(m_fluxOrders[i1]);
        ptwXYPoints *__flux = ptwXY_createFrom_Xs_Ys( nullptr, ptwXY_interpolationLinLin, ptwXY_interpolationToString( ptwXY_interpolationLinLin ), 
                10, 1e-3, 10, 10, __fluxOrder->size( ), __fluxOrder->energies( ), __fluxOrder->fluxes( ), 0 );
        if( __flux == nullptr ) throw Exception( "ptwXY_createFrom_Xs_Ys failed for __flux." );

        ptwXPoints *groupedFluxX = ptwXY_groupOneFunction( nullptr, __flux, boundaries, ptwXY_group_normType_none, nullptr );
        if( groupedFluxX == nullptr ) throw Exception( "ptwXY_groupOneFunction failed for groupedFluxX." );

        for( int i2 = 0; i2 < ptwX_length( nullptr, groupedFluxX ); ++i2 ) groupedFlux.push_back( ptwX_getPointAtIndex_Unsafely( groupedFluxX, i2 ) );

        ptwX_free( groupedFluxX );
        ptwXY_free( __flux );
    }
    ptwX_free( boundaries );

    return( ProcessedFlux( temperature( ), groupedFlux ) );
}

/* *********************************************************************************************************//**
 * Print the Flux to std::cout. Mainly for debugging.
 *
 * @param a_indent                  [in]    The std::string to print at the beginning.
 * @param a_outline                 [in]    If true, does not print the flux values.
 * @param a_valuesPerLine           [in]    The number of points (i.e., energy, flux pairs) to print per line.
 ***********************************************************************************************************/

void Flux::print( std::string const &a_indent, bool a_outline, int a_valuesPerLine ) const {

    std::cout << a_indent << "FLUX: label = '" << m_label << "': maximum order = " << ( size( ) - 1 ) << std::endl;
    if( a_outline ) return;
    for( std::vector<Flux_order>::const_iterator iter = m_fluxOrders.begin( ); iter < m_fluxOrders.end( ); ++iter )
        iter->print( a_valuesPerLine );
}


/*! \class Fluxes_from_bdfls
 * Specifies the flux data for a specified Legendre order (see class Flux).
 */

/* *********************************************************************************************************//**
 * Reads in fluxes from a *bdfls* file as a list of Flux instances.
 *
 * @param a_fileName                [in]    The bdfls file name.
 * @param a_temperature_MeV         [in]    The temperature to assign to the read fluxes.
 ***********************************************************************************************************/

Fluxes_from_bdfls::Fluxes_from_bdfls( std::string const &a_fileName, double a_temperature_MeV = 0 ) {

    initialize( a_fileName.c_str( ), a_temperature_MeV );
}

/* *********************************************************************************************************//**
 * Reads in fluxes from a *bdfls* file as a list of Flux instances.
 *
 * @param a_fileName                [in]    The bdfls file name.
 * @param a_temperature_MeV         [in]    The temperature to assign to the read fluxes.
 ***********************************************************************************************************/

Fluxes_from_bdfls::Fluxes_from_bdfls( char const *a_fileName, double a_temperature_MeV = 0 ) {

    initialize( a_fileName, a_temperature_MeV );
}

/* *********************************************************************************************************//**
 * Used by constructors to do most of the work.
 *
 * @param a_fileName                [in]    The bdfls file name.
 * @param a_temperature_MeV         [in]    The temperature to assign to the read fluxes.
 ***********************************************************************************************************/

void Fluxes_from_bdfls::initialize( char const *a_fileName, double a_temperature_MeV ) {

    char buffer[132], *pEnd, cValue[16];
    long numberOfValuesInOrders[16];
    FILE *fIn = fopen( a_fileName, "r" );
    if( fIn == nullptr ) throw Exception( "Fluxes_from_bdfls::initialize: Could not open bdfls file." );

    while( true ) {                 // Skip over groups.
        if( fgets( buffer, 132, fIn ) == nullptr ) throw Exception( "Fluxes_from_bdfls::initialize: fgets failed for fid." );
        if( strlen( buffer ) > 73 ) {
            if( buffer[72] == '1' ) break;
        }
    }

    while( true ) {
        int fid( -1 );
        if( fgets( buffer, 132, fIn ) == nullptr ) throw Exception( "Fluxes_from_bdfls::initialize: fgets failed for fid." );
        if( strlen( buffer ) > 73 ) {
            if( buffer[72] == '1' ) break;
        }
        fid = (int) strtol( buffer, &pEnd, 10 );
        if( fid == -1 ) throw Exception( "Fluxes_from_bdfls::initialize: converting fid to long failed." );
        std::string label( LLNL_fidToLabel( fid ) );
        Flux flux( label, a_temperature_MeV );

        long maximumFluxOrder( -1 );
        if( fgets( buffer, 132, fIn ) == nullptr ) throw Exception( "Fluxes_from_bdfls::initialize: fgets failed for maximumFluxOrder." );
        maximumFluxOrder = strtol( buffer, &pEnd, 10 );
        if( maximumFluxOrder == -1 ) throw Exception( "Fluxes_from_bdfls::initialize: converting maximumFluxOrder to long failed." );
        if( maximumFluxOrder >= (long) ( sizeof( numberOfValuesInOrders ) / sizeof( numberOfValuesInOrders[0] ) ) )
            throw Exception( "Fluxes_from_bdfls::initialize: need to increase size of numberOfValuesInOrders" );

        for( long order = 0; order <= maximumFluxOrder; ++order ) {
            numberOfValuesInOrders[order] = -1;
            if( fgets( buffer, 132, fIn ) == nullptr ) throw Exception( "Fluxes_from_bdfls::initialize: fgets failed for maximumFluxOrders." );
            numberOfValuesInOrders[order] = strtol( buffer, &pEnd, 10 );
            if( numberOfValuesInOrders[order] == -1 ) throw Exception( "Fluxes_from_bdfls::initialize: converting numberOfValuesInOrders[order] to long failed." );
            numberOfValuesInOrders[order] /= 2;
        }

        for( long order = 0; order <= maximumFluxOrder; ++order ) {
            long index = 0, numberOfValuesInOrder = 2 * numberOfValuesInOrders[order];
            std::vector<double> energiesAndFluxes( numberOfValuesInOrder );
            while( numberOfValuesInOrder > 0 ) {
                long i1, n1 = 6;
                if( numberOfValuesInOrder < 6 ) n1 = numberOfValuesInOrder;
                if( fgets( buffer, 132, fIn ) == nullptr ) throw Exception( "Fluxes_from_bdfls::initialize: fgets failed for energies/fluxes." );
                for( i1 = 0; i1 < n1; ++i1, ++index ) {
                    strncpy( cValue, &buffer[12*i1], 12 );
                    cValue[12] = 0;
                    energiesAndFluxes[index] = strtod( cValue, &pEnd );
                }
                numberOfValuesInOrder -= n1;
            }

            std::vector<double> energies( numberOfValuesInOrders[order] );
            std::vector<double> fluxes( numberOfValuesInOrders[order] );
            for( index = 0; index < numberOfValuesInOrders[order]; ++index ) {
                energies[index] = energiesAndFluxes[2*index];
                fluxes[index] = energiesAndFluxes[2*index+1];
            }

            Flux_order flux_order( order, energies, fluxes );
            flux.addFluxOrder( flux_order );
        }
        m_fluxes.push_back( flux );
    }

    fclose( fIn );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Fluxes_from_bdfls::~Fluxes_from_bdfls( ) {

}

/* *********************************************************************************************************//**
 * Returns the Flux instance whose *fid* is *a_fid*.
 *
 * @param a_fid             [in]    The *fid* of the Flux to return.
 * @return                          Returns the Flux whose *fid* is *a_fid*.
 ***********************************************************************************************************/

Flux Fluxes_from_bdfls::getViaFID( int a_fid ) const {

    std::string label( LLNL_fidToLabel( a_fid ) );

    for( int if1 = 0; if1 < (int) m_fluxes.size( ); ++if1 ) {
        if( m_fluxes[if1].label( ) == label ) return( m_fluxes[if1] );
    }
    throw Exception( "Fluxes_from_bdfls::getViaFID: fid not found." );
}

/* *********************************************************************************************************//**
 * Returns the 3-d function flux f(T,E,mu) whose *fid* is *a_fid*. In f(T,E,mu), T is the temperature, E is the projectile'e energy
 * and mu is the cos(theta) where theta is measured relative to the projectile velocity.
 *
 * @param a_fid             [in]    The *fid* of the Flux to return.
 * @return                          Returns the 3-d function flux f(T,E,mu).
 ***********************************************************************************************************/

Functions::XYs3d *Fluxes_from_bdfls::get3dViaFID( int a_fid ) const {

    Flux flux = getViaFID( a_fid );

    Axes axes;
    axes.append( new Axis( 3, "temperature", "MeV/k" ) );
    axes.append( new Axis( 2, "energy_in", "MeV" ) );
    axes.append( new Axis( 1, "mu", "" ) );
    axes.append( new Axis( 0, "flux", "1/s" ) );

    Functions::XYs3d *xys3d = new Functions::XYs3d( axes, ptwXY_interpolationLinLin, 0, flux.temperature( ) );
    Functions::XYs2d *xys2d = new Functions::XYs2d( axes, ptwXY_interpolationLinLin );

    xys3d->setLabel( flux.label( ) );

    std::vector<double> const &energies = flux[0].v_energies( );
    for( std::size_t i1 = 0; i1 < energies.size( ); ++i1 ) {
        Functions::Legendre1d *legendre1d = new Functions::Legendre1d( axes, 0, energies[i1] );
        std::vector<double> &coefficients = legendre1d->coefficients( );

        for( int i2 = 0; i2 < flux.size( ); ++i2 ) coefficients.push_back( flux[i2].fluxes( )[i1] );
        xys2d->append( legendre1d );
    }
    xys3d->append( xys2d );

    return( xys3d );
}

/* *********************************************************************************************************//**
 * Returns a list of *fid* strings (i.e., labels) for the Flux's read in from *bdfls* file.
 *
 * @return                          The list of *fid*'s.
 ***********************************************************************************************************/

std::vector<std::string> Fluxes_from_bdfls::labels( ) const {

    int size = (int) m_fluxes.size( );
    std::vector<std::string> _labels( size );

    for( int if1 = 0; if1 < size; ++if1 ) _labels[if1] = m_fluxes[if1].label( );
    return( _labels );
}

/* *********************************************************************************************************//**
 * Returns a list of integer *fid* for the Flux's read in from *bdfls* file.
 *
 * @return                          The list of *fid*'s.
 ***********************************************************************************************************/

std::vector<int> Fluxes_from_bdfls::FIDs( ) const {

    int size = (int) m_fluxes.size( );
    std::vector<int> fids( size );
    char *e;

    for( int if1 = 0; if1 < size; ++if1 ) {
        fids[if1] = (int) strtol( &(m_fluxes[if1].label( ).c_str( )[9]), &e, 10 );
    }
    return( fids );
}

/* *********************************************************************************************************//**
 * Print the list of Flux's to std::cout. Mainly for debugging.
 *
 * @param a_outline                 [in]    Passed to other *print* methods.
 * @param a_valuesPerLine           [in]    Passed to other *print* methods.
 ***********************************************************************************************************/

void Fluxes_from_bdfls::print( bool a_outline, int a_valuesPerLine ) const {

    int nfs = (int) m_fluxes.size( );

    std::cout << "BDFLS FLUXes: number of fluxes = " << nfs << std::endl;
    for( int if1 = 0; if1 < nfs ; ++if1 ) m_fluxes[if1].print( "  ", a_outline, a_valuesPerLine );
}

}

/* *********************************************************************************************************//**
 * Convert a flux in the form of a Function3dForm into those needed by Settings methods. Currently only works for
 * XYs3d which contains a list of XYs2d which each contain a list of Legendre instances.
 *
 * @param a_function3d          [in]    A Function3dForm instance.
 * @return                              The list of Transporting::Flux * instances.
 ***********************************************************************************************************/

std::vector<Transporting::Flux> settingsFluxesFromFunction3d( Functions::Function3dForm const &a_function3d ) {

    if( a_function3d.type( ) != FormType::XYs3d ) throw Exception( "Currently, only a 3d function of type XYs3d is supported." );

    Functions::XYs3d const &xys3d = static_cast<Functions::XYs3d const &>( a_function3d );
    std::vector<Functions::Function2dForm *> const function2ds = xys3d.function2ds( );
    std::vector<Transporting::Flux> fluxes;

    for( std::size_t i1 = 0; i1 < function2ds.size( ); ++i1 ) {
        Functions::Function2dForm const &function2d = *function2ds[i1];

        if( function2d.type( ) != FormType::XYs2d ) throw Exception( "Currently, only a 2d function of type XYs2d is supported for flux f(E,mu)." );

        Functions::XYs2d const &xys2d = static_cast<Functions::XYs2d const &>( function2d );
        std::vector<Functions::Function1dForm *> const &function1ds = xys2d.function1ds( );

        Transporting::Flux flux( a_function3d.label( ), xys2d.outerDomainValue( ) );
        std::size_t maxOrder = 0;
        std::vector<double> energies;
        std::vector< std::vector<double> > fluxMatrix;

        for( std::size_t i2 = 0; i2 < function1ds.size( ); ++i2 ) {
            Functions::Function1dForm const &function1d = *function1ds[i2];

            if( function1d.type( ) != FormType::Legendre1d ) throw Exception( "Currently, only a 1d function of type Legendre1d is supported for flux f(mu)." );

            Functions::Legendre1d const &legendre1d = static_cast<Functions::Legendre1d const &>( function1d );

            energies.push_back( legendre1d.outerDomainValue( ) );
            std::vector<double> &coefficients = const_cast< std::vector<double> &>( legendre1d.coefficients( ) );
            if( maxOrder < coefficients.size( ) ) maxOrder = coefficients.size( );
            fluxMatrix.push_back( coefficients );
        }

        for( std::size_t order = 0; order < maxOrder; ++order ) {
            std::vector<double> energyFluxAtOrder;

            for( std::size_t i2 = 0; i2 < function1ds.size( ); ++i2 ) {
                energyFluxAtOrder.push_back( 0.0 );
                if( order < fluxMatrix[i2].size( ) ) energyFluxAtOrder[i2] = fluxMatrix[i2][order];
            }

            Transporting::Flux_order fluxOrder( order, energies, energyFluxAtOrder );
            flux.addFluxOrder( fluxOrder );
        }

        fluxes.push_back( flux );
    }

    return( fluxes );
}

}
