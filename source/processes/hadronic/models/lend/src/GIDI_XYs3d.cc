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

/*! \class XYs3d
 * Class for the GNDS <**XYs3d**> node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

XYs3d::XYs3d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function3dForm( a_construction, a_node, a_setupInfo, FormType::XYs3d, a_parent ),
        m_interpolationQualifier( a_node.attribute_as_string( GIDI_interpolationQualifierChars ) ) {

    if( a_setupInfo.m_formatVersion.format( ) != GNDS_formatVersion_1_10Chars ) {
        data2dListParse( a_construction, a_node.child( GIDI_function2dsChars ), a_setupInfo, m_function2ds );
        checkOuterDomainValues2d( m_function2ds, m_Xs );
        return;                                     // Need to add uncertainty parsing.
    }

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        std::string name( child.name( ) );

        if( name == GIDI_axesChars ) continue;
        if( name == GIDI_uncertaintyChars ) continue;

        Function2dForm *_form = data2dParse( a_construction, child, a_setupInfo, nullptr );
        if( _form == nullptr ) throw Exception( "XYs3d::XYs3d: data2dParse returned nullptr." );
        if( m_Xs.size( ) > 0 ) {
            if( _form->outerDomainValue( ) <= m_Xs[m_Xs.size( )-1] ) throw Exception( "XYs3d::XYs3d: next outerDomainValue <= current outerDomainValue." );
        }
        m_Xs.push_back( _form->outerDomainValue( ) );
        m_function2ds.push_back( _form );
    }

}

/* *********************************************************************************************************//**
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 ***********************************************************************************************************/

XYs3d::XYs3d( Axes const &a_axes, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue ) :
        Function3dForm( GIDI_XYs3dChars, FormType::XYs3d, a_axes, a_interpolation, a_index, a_outerDomainValue ),
        m_interpolationQualifier( "" ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

XYs3d::~XYs3d( ) {

    for( std::vector<Function2dForm *>::iterator iter = m_function2ds.begin( ); iter < m_function2ds.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Returns the domain minimum for the instance.
 *
 * @return          The domain minimum for the instance.
 ***********************************************************************************************************/

double XYs3d::domainMin( ) const {

    if( m_Xs.size( ) == 0 ) throw Exception( "XYs3d::domainMin: XYs3d has no 2d-functions" );
    return( m_Xs[0] );
}

/* *********************************************************************************************************//**
 * Returns the domain maximum for the instance.
 *
 * @return              The domain maximum for the instance.
 ***********************************************************************************************************/

double XYs3d::domainMax( ) const {

    if( m_Xs.size( ) == 0 ) throw Exception( "XYs3d::domainMax: XYs3d has no 2d-functions" );
    return( m_Xs[m_Xs.size( )-1] );
}

/* *********************************************************************************************************//**
 * Returns the value of the function *f(x3,x2,x1)* at the specified point *a_x3*, *a_x2* and *a_x1*.
 *
 * @param a_x3              [in]    The value of the **x3** axis.
 * @param a_x2              [in]    The value of the **x2** axis.
 * @param a_x1              [in]    The value of the **x1** axis.
 * @return                          The value of the function evaluated at *a_x3*, *a_x2* and *a_x1*.
 ***********************************************************************************************************/

double XYs3d::evaluate( double a_x3, double a_x2, double a_x1 ) const {

    if( m_Xs.size( ) == 0 ) throw Exception( "XYs3d::evaluate: XYs3d has no 2d functions." );

    long iX3 = binarySearchVector( a_x3, m_Xs );

    if( iX3 < 0 ) {
        if( iX3 == -1 ) {       /* x3 > last value of Xs. */
            return( m_function2ds.back( )->evaluate( a_x2, a_x1 ) );
        }
        return( m_function2ds[0]->evaluate( a_x2, a_x1 ) ); /* x3 < first value of Xs. */
    }

// Currently does not interpolate;
    return( m_function2ds[static_cast<std::size_t>(iX3)]->evaluate( a_x2, a_x1 ) );
}

/* *********************************************************************************************************//**
 * Appends *a_function2d* to the end of *this*.
 * 
 * @param a_function2d      [in]    The 2d-function to append to *this*.
 ***********************************************************************************************************/

void XYs3d::append( Function2dForm *a_function2d ) {

    if( m_function2ds.size( ) > 0 ) {
        if( a_function2d->outerDomainValue( ) <= m_Xs.back( ) ) throw Exception( "XYs3d::append: outerDomainValue not greater than current maximum outer domain's value." );
    }

    m_Xs.push_back( a_function2d->outerDomainValue( ) );
    m_function2ds.push_back( a_function2d );
    a_function2d->setAncestor( this );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions2d container.
 ***********************************************************************************************************/

void XYs3d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const {

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

    if( m_interpolationQualifier != "" ) attributes = a_writeInfo.addAttribute( GIDI_interpolationQualifierChars, m_interpolationQualifier );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    axes( ).toXMLList( a_writeInfo, indent2 );
    for( std::vector<Function2dForm *>::const_iterator iter = m_function2ds.begin( ); iter != m_function2ds.end( ); ++iter ) (*iter)->toXMLList_func( a_writeInfo, indent2, true, false );
    a_writeInfo.addNodeEnder( moniker( ) );
}

}               // End namespace Functions.

}               // End namespace GIDI.
