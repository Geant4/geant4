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

/*! \class Gridded1d
 * Class for the GNDS <**gridded1d**> node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Gridded1d::Gridded1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::gridded1d, a_parent ) {

    Grid const *axis = dynamic_cast<Grid const *>( axes( )[0] );
    m_grid = axis->data( ).vector();

    parseFlattened1d( a_construction, a_node.child( GIDI_arrayChars ), a_setupInfo, m_data );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Gridded1d::~Gridded1d( ) {

}

/* *********************************************************************************************************//**
 * Only for internal use. Called by ProtareTNSL instance to zero the lower energy multi-group data covered by the ProtareSingle that
 * contains the TNSL data covers the lower energy multi-group data.
 *
 * @param a_maxTNSL_index           [in]    All elements up to *a_maxTNSL_index* exclusive are zero-ed.
 ***********************************************************************************************************/

void Gridded1d::modifiedMultiGroupElasticForTNSL( int a_maxTNSL_index ) {

    m_data.setToValueInFlatRange( 0, a_maxTNSL_index, 0.0 );
}

/* *********************************************************************************************************//**
 * Returns the value of the function at the point *a_x1*.
 * Currently not implemented.
 *
 * @param a_x1              [in]    The is ignored.
 * @return                          The value of the function at the point *a_x1*.
 ***********************************************************************************************************/

double Gridded1d::evaluate( LUPI_maybeUnused double a_x1 ) const {

    throw Exception( "Gridded1d::evaluate: not implement." );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions1d container.
 ***********************************************************************************************************/

void Gridded1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const {

// BRB. This is not correct as it is not converted to a flattened array.

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string indent3 = a_writeInfo.incrementalIndent( indent2 );
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

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    axes( ).toXMLList( a_writeInfo, indent2 );

    attributes = a_writeInfo.addAttribute( GIDI_shapeChars, size_t_ToString( m_data.size( ) ) );
    attributes += a_writeInfo.addAttribute( "compression", "flattened" );
    a_writeInfo.addNodeStarter( indent2, GIDI_arrayChars, attributes );

    std::vector<double> doubles;
    doubles.reserve( m_data.size( ) );
    std::size_t i1, i2;
    for( i1 = 0; i1 < m_data.size( ); ++i1 ) {
        if( m_data[i1] != 0.0 ) break;
    }
    for( i2 = m_data.size( ); i2 > i1; --i2 ) {
        if( m_data[i2-1] != 0.0 ) break;
    }
    std::size_t start( i1 );
    if( start == m_data.size( ) ) start = 0;
    a_writeInfo.push_back( indent3 + "<values valueType=\"Integer32\" label=\"starts\">" + size_t_ToString( start ) + "</values>" );
    for( ; i1 < i2; ++i1 ) doubles.push_back( m_data[i1] );
    a_writeInfo.push_back( indent3 + "<values valueType=\"Integer32\" label=\"lengths\">" + size_t_ToString( doubles.size( ) ) + "</values>" );

    doublesToXMLList( a_writeInfo, indent3, doubles );
    a_writeInfo.addNodeEnder( GIDI_arrayChars );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/* *********************************************************************************************************//**
 * This method writes *this* to a *a_file*.
 *
 * @param   a_file          [in]    The C FILE instance to write the data to.
 * @param   a_format        [in]    The format string passed to each region's write method.
 ***********************************************************************************************************/

void Gridded1d::write( FILE *a_file, std::string const &a_format ) const {

    std::size_t index = 0;
    char const *fmt = a_format.c_str( );

    for( ; index < m_data.size( ); ++index ) {
        fprintf( a_file, fmt, m_grid[index], m_data[index] );
    }
    if( m_data.size( ) > 0 ) printf( fmt, m_grid[index], m_data[index-1] );
}

}               // End namespace Functions.

}               // End namespace GIDI.
