/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <GIDI.hpp>

namespace GIDI {

namespace Table {

/*! \class Table
 * Class for the GNDS **table** node.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

Table::Table( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::table ),
        m_rows( static_cast<std::size_t>( a_node.attribute_as_int( GIDI_rowsChars ) ) ),
        m_columns( static_cast<std::size_t>( a_node.attribute_as_int( GIDI_columnsChars ) ) ),
        m_storageOrder( a_node.attribute_as_string( GIDI_storageOrderChars ) ),
        m_columnHeaders( a_construction, GIDI_columnHeadersChars, GIDI_indexChars, a_node, a_setupInfo, PoPI::Database( ), PoPI::Database( ), parseColumnHeaders, nullptr ),
        m_data( a_construction, a_node.child( GIDI_dataChars ), a_setupInfo ) {

    if( m_storageOrder == "" ) m_storageOrder = GIDI_rowMajorChars;

    m_columnHeaders.setAncestor( this );
    m_data.setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Table::~Table( ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Table::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    std::string attributes;
    attributes += a_writeInfo.addAttribute( GIDI_rowsChars, intToString( static_cast<int>( m_rows ) ) );
    attributes += a_writeInfo.addAttribute( GIDI_columnsChars, intToString( static_cast<int>( m_columns ) ) );
    if( m_storageOrder != GIDI_rowMajorChars ) attributes += a_writeInfo.addAttribute( GIDI_storageOrderChars, m_storageOrder );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    m_columnHeaders.toXMLList( a_writeInfo, indent2 );
    m_data.toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class Column
 * Class for the GNDS **column** node that is a child of the **columnHeaders** node.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Column::Column( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Form( a_node, a_setupInfo, FormType::column, a_parent ),
        m_index( a_node.attribute_as_string( GIDI_indexChars ) ),
        m_name( a_node.attribute_as_string( GIDI_nameChars ) ),
        m_unit( a_node.attribute_as_string( GIDI_unitChars ) ),
        m_types( a_node.attribute_as_string( GIDI_typesChars ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Column::~Column( ) {

}

/* *********************************************************************************************************//**
 * Set the *m_keyValue* per the *a_keyName* name. This method assumes that *a_keyName* is "label". Otherwise, it executes a throw.
 *
 * @param a_keyName             [in]    The name of the key whose value is set.
 ***********************************************************************************************************/

void Column::setKeyValue( std::string const &a_keyName ) const {

    if( a_keyName != GIDI_indexChars ) throw Exception( "Form::setKeyValue: unsupported keyname \""  + a_keyName + "\"." );

    m_keyValue = m_index;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Column::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    
    std::string attributes;
    attributes += a_writeInfo.addAttribute( GIDI_indexChars, m_index );
    attributes += a_writeInfo.addAttribute( GIDI_nameChars, m_name );
    if( m_unit != "" ) attributes += a_writeInfo.addAttribute( GIDI_unitChars, m_unit );
    if( m_types != "" ) attributes += a_writeInfo.addAttribute( GIDI_typesChars, m_types );
    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}

/*! \class Data
 * Class for the GNDS **data** node that is a child of the **table** node.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

Data::Data( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, LUPI_maybeUnused SetupInfo &a_setupInfo ) :
        GUPI::Ancestry( GIDI_dataChars ),
        m_sep( a_node.attribute_as_string( GIDI_sepChars ) ),
        m_body( a_node.text().get() ) {

    if( m_sep == "" ) m_sep = " ";
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Data::~Data( ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Data::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;
    if( m_sep != " " ) attributes += a_writeInfo.addAttribute( GIDI_sepChars, m_sep );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    a_writeInfo.push_back( m_body );

    a_writeInfo.addNodeEnder( moniker( ) );
}

}               // End namespace Table.

}               // End namespace GIDI.
