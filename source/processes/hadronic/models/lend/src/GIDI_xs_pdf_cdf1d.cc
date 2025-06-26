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

#define GIDI_xsChars "xs"
#define GIDI_pdfChars "pdf"
#define GIDI_cdfChars "cdf"

/*! \class Xs_pdf_cdf1d
 * Class for the GNDS <**xs_pdf_cdf1d**> node.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Xs_pdf_cdf1d::Xs_pdf_cdf1d( ) :
        Function1dForm( GIDI_xs_pdf_cdf1dChars, FormType::XYs1d, Axes(), ptwXY_interpolationLinLin, 0, 0.0 ) {

}

/* *********************************************************************************************************//**
 *
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_interpolation       [in]    The interpolation flag.
 * @param a_index               [in]    If imbedded in a two dimensional function, the index of this instance.
 * @param a_outerDomainValue    [in]    If imbedded in a two dimensional function, the domain value for *x2*.
 * @param a_Xs                  [in]    List of x1 values.
 * @param a_pdf                 [in]    The pdf evaluated at the x1 values.
 * @param a_cdf                 [in]    The pdf evaluated at the x1 values.
 ***********************************************************************************************************/

Xs_pdf_cdf1d::Xs_pdf_cdf1d( Axes const &a_axes, ptwXY_interpolation a_interpolation, std::vector<double> const &a_Xs, 
                std::vector<double> const &a_pdf, std::vector<double> const &a_cdf, int a_index, double a_outerDomainValue ) :
        Function1dForm( GIDI_xs_pdf_cdf1dChars, FormType::xs_pdf_cdf1d, a_axes, a_interpolation, a_index, a_outerDomainValue ),
        m_xs( a_Xs ),
        m_pdf( a_pdf ),
        m_cdf( a_cdf ) {

}

/* *********************************************************************************************************//**
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Xs_pdf_cdf1d::Xs_pdf_cdf1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::xs_pdf_cdf1d, a_parent ) {

        nf_Buffer<double> buffer;
        parseValuesOfDoubles( a_construction, a_node.child( GIDI_xsChars ).child( GIDI_valuesChars ), a_setupInfo, buffer );
        m_xs = buffer.vector();
        parseValuesOfDoubles( a_construction, a_node.child( GIDI_pdfChars ).child( GIDI_valuesChars ), a_setupInfo, buffer );
        m_pdf = buffer.vector();
        parseValuesOfDoubles( a_construction, a_node.child( GIDI_cdfChars ).child( GIDI_valuesChars ), a_setupInfo, buffer );
        m_cdf = buffer.vector();
}

/* *********************************************************************************************************//**
 * *********************************************************************************************************/

Xs_pdf_cdf1d::~Xs_pdf_cdf1d( ) {

}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the members of *this* to those of *a_rhs* except for those
 * not set by base classes.
 *
 * @param a_rhs                     [in]    Instance whose member are used to set the members of *this*.
 ***********************************************************************************************************/

Xs_pdf_cdf1d &Xs_pdf_cdf1d::operator=( Xs_pdf_cdf1d const &a_rhs ) {

    if( this != &a_rhs ) {
        Function1dForm::operator=( a_rhs );
        m_xs = a_rhs.Xs( );
        m_pdf = a_rhs.pdf( );
        m_cdf = a_rhs.cdf( );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * The value of *pdf* at the point *a_x1*.
 * Currently not implemented.
 *
 * @param a_x1          [in]    The point for the *x1* axis.
 * @return                      The value of the function at the point *a_x1*.
 ***********************************************************************************************************/

double Xs_pdf_cdf1d::evaluate( LUPI_maybeUnused double a_x1 ) const {

    throw Exception( "Xs_pdf_cdf1d::evaluate: not implemented." );
}

/* *********************************************************************************************************//**
 * This methods returns an XYs1d representation of the pdf of *this*. The calling function owns the created instance and is responible
 * for freeing it.
 *
 * @param   a_asLinlin          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 * @param   a_accuracy          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 * @param   a_lowerEps          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 * @param   a_upperEps          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 *
 * @return                                  A pointer to an  XYs1d instance that must be freed by the calling function.
 ***********************************************************************************************************/

XYs1d *Xs_pdf_cdf1d::asXYs1d( LUPI_maybeUnused bool a_asLinlin, LUPI_maybeUnused double a_accuracy, LUPI_maybeUnused double a_lowerEps, LUPI_maybeUnused double a_upperEps ) const {

    return( new XYs1d( axes( ), ptwXY_interpolationLinLin, m_xs, m_pdf ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions1d container.
 ***********************************************************************************************************/

void Xs_pdf_cdf1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const {

    std::string attributes;

    if( a_embedded ) {
        attributes += a_writeInfo.addAttribute( GIDI_outerDomainValueChars, LUPI::Misc::doubleToShortestString( outerDomainValue( ) ) ); }
    else {
        if( a_inRegions ) {
            attributes = a_writeInfo.addAttribute( GIDI_indexChars, intToString( index( ) ) ); }
        else {
            if( label( ) != "" ) attributes = a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
        }
    }

    if( interpolation( ) != ptwXY_interpolationLinLin ) attributes += a_writeInfo.addAttribute( GIDI_interpolationChars, interpolationString( ) );

    std::string xml = a_writeInfo.nodeStarter( a_indent, moniker( ), attributes );
    xml += nodeWithValuesToDoubles( a_writeInfo,  GIDI_xsChars, m_xs );
    xml += nodeWithValuesToDoubles( a_writeInfo, GIDI_pdfChars, m_pdf );
    xml += nodeWithValuesToDoubles( a_writeInfo, GIDI_cdfChars, m_cdf );
    xml += a_writeInfo.nodeEnder( moniker( ) );

    a_writeInfo.push_back( xml );
}

}               // End namespace Functions.

}               // End namespace GIDI.
