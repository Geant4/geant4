/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

namespace Functions {

/*! \class Constant1d
 * Class for the GNDS <**constant1d**> node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_value               [in]     The GNDS **value** for *this*.
 * @param a_domainMin           [in]     The minimum value for the domain.
 * @param a_domainMax           [in]     The maximum value for the domain.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 ***********************************************************************************************************/

Constant1d::Constant1d( Axes const &a_axes, double a_value, double a_domainMin, double a_domainMax, int a_index, double a_outerDomainValue ) :
        Function1dForm( GIDI_constant1dChars, FormType::constant1d, a_axes, ptwXY_interpolationLinLin, a_index, a_outerDomainValue ),
        m_value( a_value ),
        m_domainMin( a_domainMin ),
        m_domainMax( a_domainMax ) {

}

/* *********************************************************************************************************//**
 *
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Constant1d::Constant1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::constant1d, a_parent ),
        m_value( a_node.attribute( GIDI_valueChars ).as_double( ) ),
        m_domainMin( a_node.attribute( GIDI_domainMinChars ).as_double( ) ),
        m_domainMax( a_node.attribute( GIDI_domainMaxChars ).as_double( ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Constant1d::~Constant1d( ) {

}

/* *********************************************************************************************************//**
 * Returns the value of the constant function.
 *
 * @param a_x1              [in]    This is ignored a the function is a constant.
 * @return                          The value of the constant.
 ***********************************************************************************************************/

double Constant1d::evaluate( LUPI_maybeUnused double a_x1 ) const {

// FIXME - Do we need to check domain?
    return( m_value );
}

/* *********************************************************************************************************//**
 * Evaluates *this* at the X-values in *a_Xs*[*a_offset*:] and adds the results to *a_results*[*a_offset*:].
 * *a_Xs* and *a_results* must be the same size otherwise a throw is executed.
 *
 * @param a_offset          [in]    The offset in *a_Xs* to start.
 * @param a_Xs              [in]    The list of domain values to evaluate *this* at.
 * @param a_results         [in]    The list whose values are added to by the Y-values of *this*.
 * @param a_scaleFactor     [in]    A factor applied to each evaluation before it is added to *a_results*.
 ***********************************************************************************************************/

void Constant1d::mapToXsAndAdd( std::size_t a_offset, std::vector<double> const &a_Xs, std::vector<double> &a_results, double a_scaleFactor ) const {

    if( a_Xs.size( ) != a_results.size( ) ) throw Exception( "Constant1d::mapToXsAndAdd: a_Xs.size( ) != a_results.size( )." );

    for( std::size_t index = a_offset; index < a_Xs.size( ); ++index ) {
        a_results[index] += a_scaleFactor * m_value;
    }
}

/* *********************************************************************************************************//**
 * This methods returns an XYs1d representation of *this*. The calling function owns the created instance and is responible
 * for freeing it.
 *
 * @param   a_asLinlin          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 * @param   a_accuracy          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 * @param   a_lowerEps          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 * @param   a_upperEps          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 *
 * @return                                  A pointer to an  XYs1d instance that must be freed by the calling function.
 ***********************************************************************************************************/

XYs1d *Constant1d::asXYs1d( LUPI_maybeUnused bool a_asLinlin, LUPI_maybeUnused double a_accuracy, LUPI_maybeUnused double a_lowerEps, LUPI_maybeUnused double a_upperEps ) const {

    std::vector<double> xys( 4 );

    xys[0] = m_domainMin;
    xys[1] = m_value;
    xys[2] = m_domainMax;
    xys[3] = m_value;

    return( new XYs1d( axes( ), ptwXY_interpolationLinLin, xys ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions1d container.
 ***********************************************************************************************************/

void Constant1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes;
    
    if( a_embedded ) {
        attributes += a_writeInfo.addAttribute( GIDI_outerDomainValueChars, LUPI::Misc::doubleToShortestString( outerDomainValue( ) ) ); }
    else {
        if( a_inRegions ) {
            attributes = a_writeInfo.addAttribute( GIDI_indexChars, intToString( index( ) ) ); }
        else {
            attributes = a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
        }
    }

    attributes += a_writeInfo.addAttribute( GIDI_valueChars, LUPI::Misc::doubleToShortestString( m_value ) );
    attributes += a_writeInfo.addAttribute( GIDI_domainMinChars, LUPI::Misc::doubleToShortestString( m_domainMin ) );
    attributes += a_writeInfo.addAttribute( GIDI_domainMaxChars, LUPI::Misc::doubleToShortestString( m_domainMax ) );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    axes( ).toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

}               // End namespace Functions.

}               // End namespace GIDI.
