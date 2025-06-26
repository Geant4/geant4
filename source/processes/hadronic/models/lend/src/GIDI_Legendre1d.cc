/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>

#include <nf_Legendre.h>
#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

namespace Functions {

/*! \class Legendre1d
 * Class for the GNDS <**Legendre**> node.
 */

/* *********************************************************************************************************//**
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 ***********************************************************************************************************/

Legendre1d::Legendre1d( Axes const &a_axes, int a_index, double a_outerDomainValue ) :
        Function1dForm( GIDI_LegendreChars, FormType::Legendre1d, a_axes, ptwXY_interpolationLinLin, a_index, a_outerDomainValue ) {

}

/* *********************************************************************************************************//**
 *
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Legendre1d::Legendre1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::Legendre1d, a_parent ) {

    nf_Buffer<double> coeff;
    parseValuesOfDoubles( a_construction, a_node.child( GIDI_valuesChars ), a_setupInfo, coeff );
    m_coefficients = coeff.vector();
}

/* *********************************************************************************************************//**
 * The Legendre1d copy constructor.
 *
 * @param a_Legendre1d          [in]    The Legendre1d instance to copy.
 ***********************************************************************************************************/

Legendre1d::Legendre1d( Legendre1d const &a_Legendre1d ) :
        Function1dForm( GIDI_LegendreChars, FormType::Legendre1d, a_Legendre1d.axes( ), ptwXY_interpolationLinLin, 0, 0.0 ),
        m_coefficients( a_Legendre1d.coefficients( ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Legendre1d::~Legendre1d( ) {

}

/* *********************************************************************************************************//**
 * Returns the value of the function evaluated at the specified projectile's energy.
 * Currently not implemented.
 *
 * @param a_x1              [in]    The projectile's energy.
 * @return                          The value of the function evaluated at *a_x1*.
 ***********************************************************************************************************/

double Legendre1d::evaluate( LUPI_maybeUnused double a_x1 ) const {

    throw Exception( "Legendre1d::evaluate: not implemented." );
}

/* *********************************************************************************************************//**
 * This methods returns an XYs1d representation of *this*. The calling function owns the created instance and is responible
 * for freeing it.
 *
 * @param   a_asLinlin          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 * @param   a_accuracy          [in]    The accuracy use to convert the data to lin=lin interpolation if needed. This argument is not needed or used for this cl
 * @param   a_lowerEps          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 * @param   a_upperEps          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 *
 * @return                                  A pointer to an  XYs1d instance that must be freed by the calling function.
 ***********************************************************************************************************/

XYs1d *Legendre1d::asXYs1d( LUPI_maybeUnused bool a_asLinlin, double a_accuracy, LUPI_maybeUnused double a_lowerEps, LUPI_maybeUnused double a_upperEps ) const {

    int size1 = static_cast<int>( m_coefficients.size( ) );
    nf_Legendre *legendre1 = nf_Legendre_new( nullptr, 0, size1 - 1, const_cast<double *>( m_coefficients.data( ) ) );
    if( legendre1 == nullptr ) return( nullptr );

    ptwXYPoints *xys = nf_Legendre_to_ptwXY( nullptr, legendre1, a_accuracy, 12, 1 );
    nf_Legendre_free( legendre1 );
    if( xys == nullptr ) return( nullptr );

    return( new XYs1d( axes( ), xys ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions1d container.
 ***********************************************************************************************************/

void Legendre1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes;

    if( a_embedded ) {
        attributes += a_writeInfo.addAttribute( GIDI_outerDomainValueChars, LUPI::Misc::doubleToShortestString( outerDomainValue( ) ) ); }
    else {
        if( a_inRegions ) {
            attributes = a_writeInfo.addAttribute( GIDI_indexChars, intToString( index( ) ) ); }
        else {
            if( keyValue( ) != "" ) attributes = a_writeInfo.addAttribute( keyName( ), keyValue( ) );
        }
    }

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    if( !a_embedded ) axes( ).toXMLList( a_writeInfo, indent2 );

    doublesToXMLList( a_writeInfo, indent2, m_coefficients, 0, true );
    a_writeInfo.addNodeEnder( moniker( ) );
}

}               // End namespace Functions.

}               // End namespace GIDI.
