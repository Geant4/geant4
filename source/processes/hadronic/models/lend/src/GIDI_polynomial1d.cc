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

/* *********************************************************************************************************//**
 * This is the function callback used by asXYs1d to evalaute *this* at a domain point. This function is for internal use.
 *
 * @param       a_smr               [in/out]    
 * @param       a_xValue            [in]        The x-value to evaluate the polynomial at.
 * @param       a_yValue            [in]        A pointer to a double that will contained the polynomial evaluated at *a_xValue*.
 * @param       a_argList           [in]        A pointer to a list of additional arguments needed.
 *
 * @return                                      A nfu_status value.
 ***********************************************************************************************************/

static nfu_status asXYs1d_callback( LUPI_maybeUnused statusMessageReporting *a_smr, double a_xValue, double *a_yValue, void *a_argList ) {

    GIDI::Functions::Polynomial1d *polynomial1d = static_cast<GIDI::Functions::Polynomial1d *>( a_argList );

    *a_yValue = polynomial1d->evaluate( a_xValue );;
    return( nfu_Okay );
}

namespace GIDI {

namespace Functions {


/*! \class Polynomial1d
 * Class for the GNDS <**polynomial1d**> node.
 */

/* *********************************************************************************************************//**
 * @param a_axes            [in]    The axes to copy for *this*. 
 * @param a_domainMin       [in]    The minimum value for the domain.
 * @param a_domainMax       [in]    The maximum value for the domain.
 * @param a_coefficients    [in]    The coefficients representing the polynomial.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 ***********************************************************************************************************/

Polynomial1d::Polynomial1d( Axes const &a_axes, double a_domainMin, double a_domainMax, std::vector<double> const &a_coefficients, int a_index, double a_outerDomainValue ) :
        Function1dForm( GIDI_polynomial1dChars, FormType::polynomial1d, a_axes, ptwXY_interpolationLinLin, a_index, a_outerDomainValue ),
        m_domainMin( a_domainMin ),
        m_domainMax( a_domainMax ),
        m_coefficients( a_coefficients ) {

}

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Polynomial1d::Polynomial1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::polynomial1d, a_parent ),
        m_domainMin( a_node.attribute( GIDI_domainMinChars ).as_double( ) ),
        m_domainMax( a_node.attribute( GIDI_domainMaxChars ).as_double( ) ) {

  nf_Buffer<double> coeff;
    parseValuesOfDoubles( a_construction, a_node.child( GIDI_valuesChars ), a_setupInfo, coeff );
    m_coefficients = coeff.vector();
}

/* *********************************************************************************************************//**
 * The Polynomial1d copy constructor.
 *
 * @param a_polynomial1d        [in]    The Polynomial1d instance to copy.
 ***********************************************************************************************************/

Polynomial1d::Polynomial1d( Polynomial1d const &a_polynomial1d ) :
        Function1dForm( a_polynomial1d ),
        m_domainMin( a_polynomial1d.domainMin( ) ),
        m_domainMax( a_polynomial1d.domainMax( ) ),
        m_coefficients( a_polynomial1d.coefficients( ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Polynomial1d::~Polynomial1d( ) {

}

/* *********************************************************************************************************//**
 * The value of the polynomial at the point *a_x1*.
 *
 * @param a_x1          [in]    Domain value to evaluate this at.
 * @return                      The value of the polynomial at the point **a_x1**.
 ***********************************************************************************************************/

double Polynomial1d::evaluate( double a_x1 ) const {

    double _value = 0;

    if( a_x1 < m_domainMin ) return( 0.0 );
    if( a_x1 > m_domainMax ) return( 0.0 );

    for( std::vector<double>::const_reverse_iterator riter = m_coefficients.rbegin( ); riter != m_coefficients.rend( ); ++riter ) {
        _value = *riter + _value * a_x1;
    }
    return( _value );
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

void Polynomial1d::mapToXsAndAdd( std::size_t a_offset, std::vector<double> const &a_Xs, std::vector<double> &a_results, double a_scaleFactor ) const {

    if( a_Xs.size( ) != a_results.size( ) ) throw Exception( "Constant1d::mapToXsAndAdd: a_Xs.size( ) != a_results.size( )." );

    std::size_t index = 0;
    auto XsIter = a_Xs.begin( );

    for( ; XsIter != a_Xs.end( ); ++XsIter, ++index ) {
        if( a_offset == index ) break;
    }

    for( ; XsIter != a_Xs.end( ); ++XsIter, ++index ) {
        if( *XsIter >= m_domainMin ) break;
    }

    for( ; XsIter != a_Xs.end( ); ++XsIter, ++index ) {
        if( *XsIter > m_domainMax ) break;
        a_results[index] += a_scaleFactor * evaluate( *XsIter );
    }
}

/* *********************************************************************************************************//**
 * This methods returns an XYs1d representation of *this*. The calling function owns the created instance and is responible
 * for freeing it.
 *
 * @param   a_asLinlin          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 * @param   a_accuracy          [in]    The accuracy use to convert the data to lin=lin interpolation if needed.
 * @param   a_lowerEps          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 * @param   a_upperEps          [in]    This argument is not used but retained to make the methods API at same as other asXYs1d functions.
 *
 * @return                              A pointer to an  XYs1d instance that must be freed by the calling function.
 ***********************************************************************************************************/

XYs1d *Polynomial1d::asXYs1d( LUPI_maybeUnused bool a_asLinlin, double a_accuracy, LUPI_maybeUnused double a_lowerEps, LUPI_maybeUnused double a_upperEps ) const {

    XYs1d *xys1d = nullptr;

    if( m_coefficients.size( ) < 3 ) {
        double offset = 0.0, slope = 0.0;
        if( m_coefficients.size( ) > 0 ) {
            offset = m_coefficients[0];
            if( m_coefficients.size( ) == 2 ) slope = m_coefficients[1];
        }

        std::vector<double> xs( 2 );
        xs[0] = domainMin( );
        xs[1] = domainMax( );

        std::vector<double> ys( 2 );
        ys[0] = slope * xs[0] + offset;
        ys[1] = slope * xs[1] + offset;

        xys1d = new XYs1d( axes( ), ptwXY_interpolationLinLin, xs, ys ); }
    else {
        double xs[2] = { domainMin( ), domainMax( ) };
        ptwXYPoints *ptwXYPoints1 = ptwXY_createFromFunction( nullptr, 2, xs, asXYs1d_callback, const_cast<Polynomial1d *>( this ), a_accuracy, 1, 12 );
        if( ptwXYPoints1 != nullptr ) xys1d = new XYs1d( axes( ), ptwXYPoints1 );
    }

    return( xys1d );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions1d container.
 ***********************************************************************************************************/

void Polynomial1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
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

    attributes = a_writeInfo.addAttribute( GIDI_domainMinChars, LUPI::Misc::doubleToShortestString( domainMin( ) ) );
    attributes += a_writeInfo.addAttribute( GIDI_domainMaxChars, LUPI::Misc::doubleToShortestString( domainMax( ) ) );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    axes( ).toXMLList( a_writeInfo, indent2 );
    doublesToXMLList( a_writeInfo, indent2, m_coefficients );
    a_writeInfo.addNodeEnder( moniker( ) );
}

}               // End namespace Functions.

}               // End namespace GIDI.
