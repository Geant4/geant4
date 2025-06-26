/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <GUPI.hpp>

namespace GUPI {

/*! \class Suite
 * This class is used to store a list (i.e., suite) of similar type **GNDS** nodes.
*/

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Suite::Suite( std::string const &a_keyName ) :
        Ancestry( "" ),
        m_keyName( a_keyName ) {

}

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The **GNDS** moniker for the Suite instance.
 * @param a_keyName             [in]    The name of the key for elements of *this*.
 ***********************************************************************************************************/

Suite::Suite( std::string const &a_moniker, std::string const &a_keyName ) :
        Ancestry( a_moniker ),
        m_keyName( a_keyName ) {

}

/* *********************************************************************************************************//**
 * @param a_node                [in]    The HAPI::Node to be parsed and used to construct the Product.
 * @param a_keyName             [in]    The name of the key for referencing up child nodes.
 * @param a_parseSuite          [in]    This function to call to parse each sub-node.
 ***********************************************************************************************************/

Suite::Suite( HAPI::Node const &a_node, std::string const &a_keyName, GUPI_parseSuite a_parseSuite ) :
        Ancestry( a_node.name( ) ),
        m_keyName( a_keyName ) {

    parse( a_node, a_parseSuite );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Suite::~Suite( ) {

    for( std::vector<Entry *>::const_iterator iter = m_entries.begin( ); iter < m_entries.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * This methods parses all the child nodes of *a_node*.
 *
 * @param a_node                [in]    The HAPI::Node to be parsed and used to construct the Product.
 * @param a_parseSuite          [in]    This function to call to parse each sub-node.
 ***********************************************************************************************************/

void Suite::parse( HAPI::Node const &a_node, GUPI_parseSuite a_parseSuite ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        Entry *form = a_parseSuite( this, child );
        if( form != nullptr ) add( form );
    }
}

/* *********************************************************************************************************//**
 * Returns the index of the node in *this* that has keyValue *a_keyValue*.

 * @return                      [in]    The index of the node with keyValue *a_keyValue* in *this*.
 ***********************************************************************************************************/

int Suite::operator[]( std::string const &a_keyValue ) const {

    std::map<std::string, int>::const_iterator iter = m_map.find( a_keyValue );
    if( iter == m_map.end( ) ) {
        throw LUPI::Exception( "form '" + a_keyValue + "' not in database." );
    }

    return( iter->second );
}

/* *********************************************************************************************************//**
 * Adds the node *a_form* to *this*.
 *
 * @param a_form                [in]    The form to add.
 ***********************************************************************************************************/

void Suite::add( Entry *a_form ) {

    int i1 = 0;

    for( Suite::iterator iter = m_entries.begin( ); iter != m_entries.end( ); ++iter, ++i1 ) {
        if( (*iter)->keyValue( ) == a_form->keyValue( ) ) {
            m_entries[i1] = a_form;
            a_form->setAncestor( this );
            return;
        }
    }
    m_map[a_form->keyValue( )] = (int) m_entries.size( );
    m_entries.push_back( a_form );
    a_form->setAncestor( this );
}

/* *********************************************************************************************************//**
 * Returns the iterator to the node with keyValue *a_keyValue*.
 *
 * @param a_keyValue                        [in]    The keyValue of the node to find.
 *
 * @return                              The iterator to the node with keyValue *a_keyValue*.
 ***********************************************************************************************************/

Suite::iterator Suite::find( std::string const &a_keyValue ) {

    for( Suite::iterator iter = m_entries.begin( ); iter != m_entries.end( ); ++iter ) {
        if( (*iter)->keyName( ) == a_keyValue ) return( iter );
    }

    return( m_entries.end( ) );
}

/* *********************************************************************************************************//**
 * Returns the iterator to the node with keyValue *a_keyValue*.
 *
 * @param a_keyValue                        [in]    The keyValue of the node to find.
 *
 * @return                              The iterator to the node with keyValue *a_keyValue*.
 ***********************************************************************************************************/

Suite::const_iterator Suite::find( std::string const &a_keyValue ) const {

    for( Suite::const_iterator iter = m_entries.begin( ); iter != m_entries.end( ); ++iter ) {
        if( (*iter)->keyValue( ) == a_keyValue ) return( iter );
    }

    return( m_entries.end( ) );
}

/* *********************************************************************************************************//**
 * Returns a list of iterators to the nodes in *this* that have **GNDS** moniker *a_moniker*.
 *
 * @param a_moniker             [in]    The moniker to search for.
 *
 * @return                              List of iterators to the nodes in *this* that have moniker *a_moniker*.
 ***********************************************************************************************************/

std::vector<Suite::iterator> Suite::findAllOfMoniker( std::string const &a_moniker ) {

    std::vector<Suite::iterator> iters;

    for( Suite::iterator iter = m_entries.begin( ); iter != m_entries.end( ); ++iter ) {
        if( (*iter)->moniker( ) == a_moniker ) iters.push_back( iter );
    }

    return( iters );
}

/* *********************************************************************************************************//**
 * Returns a list of iterators to the nodes in *this* that have **GNDS** moniker *a_moniker*.
 *
 * @param a_moniker             [in]    The moniker to search for.
 *
 * @return                              List of iterators to the nodes in *this* that have moniker *a_moniker*.
 ***********************************************************************************************************/

std::vector<Suite::const_iterator> Suite::findAllOfMoniker( std::string const &a_moniker ) const {

    std::vector<Suite::const_iterator> iters;

    for( Suite::const_iterator iter = m_entries.begin( ); iter != m_entries.end( ); ++iter ) {
        if( (*iter)->moniker( ) == a_moniker ) iters.push_back( iter );
    }

    return( iters );
}

/* *********************************************************************************************************//**
 * Used by Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

Ancestry *Suite::findInAncestry3( std::string const &a_item ) {

    std::size_t index( a_item.find( '=' ) ), lastQuote = a_item.size( ) - 2;

    if( index == std::string::npos ) return( nullptr );
    ++index;
    if( index > lastQuote ) throw LUPI::Exception( "Suite::findInAncestry3: invalide xlink" );
    if( a_item[index] != '\'' ) throw LUPI::Exception( "Suite::findInAncestry3: invalid xlink, missing '." );
    ++index;
    if( a_item[lastQuote]  != '\'' ) throw LUPI::Exception( "Suite::findInAncestry3: invalid xlink, missing endl '." );

    std::string keyValue( a_item.substr( index, lastQuote - index ) );

    return( get<Ancestry>( keyValue ) );
}

/* *********************************************************************************************************//**
 * Used by Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

Ancestry const *Suite::findInAncestry3( std::string const &a_item ) const {

    std::size_t index( a_item.find( '=' ) ), lastQuote = a_item.size( ) - 2;

    if( index == std::string::npos ) return( nullptr );
    ++index;
    if( index > lastQuote ) throw LUPI::Exception( "Suite::findInAncestry3: invalide xlink" );
    if( a_item[index] != '\'' ) throw LUPI::Exception( "Suite::findInAncestry3: invalid xlink, missing '." );
    ++index;
    if( a_item[lastQuote]  != '\'' ) throw LUPI::Exception( "Suite::findInAncestry3: invalid xlink, missing endl '." );

    std::string keyValue( a_item.substr( index, lastQuote - index ) );

    return( get<Ancestry>( keyValue ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Suite::toXMLList( WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    if( size( ) == 0 ) return;

    std::string XMLLine( a_indent + "<" + moniker( ) + ">" );
    a_writeInfo.push_back( XMLLine );

    for( Suite::const_iterator iter = m_entries.begin( ); iter != m_entries.end( ); ++iter ) (*iter)->toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeEnder( moniker( ) );
}
 
/* *********************************************************************************************************//**
 * Prints the list of node keyValues to std::cout.
 *
 * @param a_header              [in]    A string printed before the list of keyValues is printed.
 ***********************************************************************************************************/

void Suite::printEntryLabels( std::string const &a_header ) const {

    std::cout << a_header << ": size = " << size( ) << std::endl;
    
    for( Suite::const_iterator iter = m_entries.begin( ); iter != m_entries.end( ); ++iter ) 
            std::cout << "    " << (*iter)->keyValue( ) << std::endl;
}

}
