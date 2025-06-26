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

/*! \class Grid
 * Class to store a **GNDS grid** node. 
 */

/* *********************************************************************************************************//**
 *
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the Grid.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_useSystem_strtod    [in]    Flag passed to the function nfu_stringToListOfDoubles.
 ***********************************************************************************************************/

Grid::Grid( HAPI::Node const &a_node, SetupInfo &a_setupInfo, int a_useSystem_strtod ) :
        Axis( a_node, a_setupInfo, FormType::grid ),
        m_style( a_node.attribute_as_string( GIDI_styleChars ) ),
        m_keyName( GIDI_indexChars ),
        m_keyValue( a_node.attribute_as_string( GIDI_indexChars ) ),
        m_interpolation( a_node.attribute_as_string( GIDI_interpolationChars ) ) {

    if( href( ) == "" ) {
        HAPI::Node values = a_node.first_child( );
        if( values.name( ) != std::string( GIDI_valuesChars ) ) throw Exception( "grid's first child not values" );

        m_valueType = values.attribute_as_string( GIDI_valueTypeChars );

        parseValuesOfDoubles( values, a_setupInfo, m_values, a_useSystem_strtod );
    }
}

/* *********************************************************************************************************//**
 * Copy constructor for Grid.
 *
 * @param a_grid                [in]    The Grid instance to copy.
 ***********************************************************************************************************/

Grid::Grid( Grid const &a_grid ) :
        Axis( a_grid ),
        m_style( a_grid.style( ) ),
        m_keyName( a_grid.keyName( ) ),
        m_keyValue( a_grid.keyValue( ) ),
        m_valueType( a_grid.valueType( ) ),
        m_values( a_grid.values( ) ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Grid::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes = a_writeInfo.addAttribute( GIDI_indexChars, intToString( index( ) ) );

    if( href( ) == "" ) {
        attributes += a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
        attributes += a_writeInfo.addAttribute( GIDI_unitChars, unit( ) );
        attributes += a_writeInfo.addAttribute( GIDI_styleChars, style( ) );
        a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
        doublesToXMLList( a_writeInfo, indent2, m_values.vector(), 0, true, m_valueType );
        a_writeInfo.addNodeEnder( moniker( ) ); }
    else {
        attributes += a_writeInfo.addAttribute( GIDI_hrefChars, href( ) );
        a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
    }
}

}
