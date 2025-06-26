/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <fstream>

#include "GUPI.hpp"

namespace GUPI {

/*! \class Ancestry
 * This is a base class inherit by most other classes. It allows one to construct a node's *xlink* or get another
 * node from its *xlink*.
 */

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The **GNDS** node's name (i.e., moniker).
 * @param a_attribute           [in]    Currently not used.
 ***********************************************************************************************************/

Ancestry::Ancestry( std::string const &a_moniker, std::string const &a_attribute ) :
        m_moniker( a_moniker ),
        m_ancestor( nullptr ),
        m_attribute( a_attribute ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Ancestry::~Ancestry( ) {

}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the member's of *this* to those of *a_ancestry* except for
 * the member *m_ancestor* which is set to **nullptr**.
 *
 * @param a_ancestry            [in]    Instance whose member are used to set the members of *this*.
 ***********************************************************************************************************/

Ancestry &Ancestry::operator=( Ancestry const &a_ancestry ) {

    if( this != &a_ancestry ) {
        m_moniker = a_ancestry.moniker( );
        m_ancestor = nullptr;
        m_attribute = a_ancestry.attribute( );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns the root node, ascending all parent nodes until one is found without an ancester. That node is returned.
 *
 * @return                              Returns the root node (i.e., the top level node).
 ***********************************************************************************************************/

Ancestry *Ancestry::root( ) {

    Ancestry *_root = this;

    while( _root->m_ancestor != nullptr ) _root = _root->m_ancestor;
    return( _root );
}

/* *********************************************************************************************************//**
 * Returns the root node, ascending all parent nodes until one is found without an ancester. That node is returned.
 *
 * @return                              Returns the root node (i.e., the top level node).
 ***********************************************************************************************************/

Ancestry const *Ancestry::root( ) const {

    Ancestry const *_root = this;

    while( _root->m_ancestor != nullptr ) _root = _root->m_ancestor;
    return( _root );
}

/* *********************************************************************************************************//**
 * Returns a pointer to the node whose *xlink* (i.e., *a_href*) is *a_href*.
 *
 * @param a_href                [in]    The *xlink* whose node is to be returned.
 * @return                              Returns the root node (i.e., the top level node).
 ***********************************************************************************************************/

Ancestry *Ancestry::findInAncestry( std::string const &a_href ) {

    std::vector<std::string> segments = LUPI::Misc::splitXLinkString( a_href );

    return( findInAncestry2( 0, segments ) );
}

/* *********************************************************************************************************//**
 * Returns a pointer to the node whose *xlink* (i.e., *a_href*) is *a_href*.
 *
 * @param a_href                [in]    The *xlink* whose node is to be returned.
 * @return                              Returns the root node (i.e., the top level node).
 ***********************************************************************************************************/

Ancestry const *Ancestry::findInAncestry( std::string const &a_href ) const {

    std::vector<std::string> segments = LUPI::Misc::splitXLinkString( a_href );

    return( findInAncestry2( 0, segments ) );
}

/* *********************************************************************************************************//**
 * Returns a pointer to the node whose *xlink* is defined by the *a_segments* argument. The *a_segments* is the *xlink*
 * divided into segments separated by the '/' character.
 *
 * @param a_index               [in]    An index into the *a_segments* whose segment is to be found at this level.
 * @param a_segments            [in]    The list of *xlink* segments.
 * @return                              Returns the root node (i.e., the top level node).
 ***********************************************************************************************************/

Ancestry *Ancestry::findInAncestry2( std::size_t a_index, std::vector<std::string> const &a_segments ) {

    Ancestry *item = this;

    if( a_index == a_segments.size( ) ) return( item );

    std::string segment( a_segments[a_index] );

    if( segment == "" ) {
        item = this->root( );
        ++a_index;
        if( a_segments[a_index] != item->moniker( ) ) return( nullptr ); }
    else if( segment == "." ) {
        }
    else if( segment == ".." ) {
        item = this->ancestor( ); }
    else {
        item = this->findInAncestry3( segment );
    }

    if( item == nullptr ) return( item );

    ++a_index;
    return( item->findInAncestry2( a_index, a_segments ) );
}

/* *********************************************************************************************************//**
 * Returns a pointer to the node whose *xlink* is defined by the *a_segments* argument. The *a_segments* is the *xlink*
 * divided into segments separated by the '/' character.
 *
 * @param a_index               [in]    An index into the *a_segments* whose segment is to be found at this level.
 * @param a_segments            [in]    The list of *xlink* segments.
 * @return                              Returns the root node (i.e., the top level node).
 ***********************************************************************************************************/

Ancestry const *Ancestry::findInAncestry2( std::size_t a_index, std::vector<std::string> const &a_segments ) const {

    Ancestry const *item = this;

    if( a_index == a_segments.size( ) ) return( item );

    std::string segment( a_segments[a_index] );

    if( segment == "" ) {
        item = this->root( );
        ++a_index;
        if( a_segments[a_index] != item->moniker( ) ) return( nullptr ); }
    else if( segment == "." ) {
        }
    else if( segment == ".." ) {
        item = this->ancestor( ); }
    else {
        item = this->findInAncestry3( segment );
    }

    if( item == nullptr ) return( item );

    ++a_index;
    return( item->findInAncestry2( a_index, a_segments ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *              
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/
        
LUPI_HOST void Ancestry::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_STD_STRING( m_moniker, a_buffer, a_mode );
    DATA_MEMBER_STD_STRING( m_attribute, a_buffer, a_mode );
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_ancestor = nullptr;
}

/* *********************************************************************************************************//**
 * Constructs and returns the *xlink* for *this*.
 *
 * @return          The constructed *xlink*.
 ***********************************************************************************************************/

std::string Ancestry::toXLink( ) const {

    std::string xlink( "/" + m_moniker + xlinkItemKey( ) );

    if( isRoot( ) ) return( xlink );
    return( m_ancestor->toXLink( ) + xlink );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Ancestry::toXMLList( LUPI_maybeUnused WriteInfo &a_writeInfo, LUPI_maybeUnused std::string const &a_indent ) const {

    std::cout << "Node '" << moniker( ) << "' needs toXMLList methods." << std::endl;
}

/* *********************************************************************************************************//**
 * Calls **toXMLList** and then writes the XML lines to the file "test.xml".
 ***********************************************************************************************************/

void Ancestry::printXML( ) const {

    WriteInfo writeInfo;

    toXMLList( writeInfo, "" );

    std::ofstream fileio;
    fileio.open( "test.xml" );
    for( std::list<std::string>::iterator iter = writeInfo.m_lines.begin( ); iter != writeInfo.m_lines.end( ); ++iter ) {
        fileio << *iter << std::endl;
    }
    fileio.close( );
}

/* *********************************************************************************************************//**
 * @param   a_incrementalIndent     [in]    The incremental amount of indentation a node adds to a sub-nodes indentation.
 * @param   a_valuesPerLine         [in]    The maximum number of integer or float values that are written per line before a new line is created.
 * @param   a_sep                   [in]    The separation character to use between integer and float values in a list.
 ***********************************************************************************************************/

WriteInfo::WriteInfo( std::string const &a_incrementalIndent, int a_valuesPerLine, std::string const &a_sep ) :
        m_incrementalIndent( a_incrementalIndent ),
        m_valuesPerLine( a_valuesPerLine ),
        m_sep( a_sep ) {

}

/* *********************************************************************************************************//**
 * Prints to contents the *this* to std::cout.
 ***********************************************************************************************************/

void WriteInfo::print( ) {

    for( auto line = m_lines.begin( ); line != m_lines.end( ); ++line ) std::cout << *line << std::endl;
}

}
