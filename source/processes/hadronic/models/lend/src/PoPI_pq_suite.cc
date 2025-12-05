/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "PoPI.hpp"

namespace PoPI {

/*! \class PQ_suite
 * Suite for storing the values in an physical quantity.
 */

/* *********************************************************************************************************//**
 * @param a_node                        [in]    The **HAPI::Node** node to be parsed.
 ***********************************************************************************************************/

PQ_suite::PQ_suite( HAPI::Node const &a_node ) :
        m_label( a_node.name( ) ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        std::string name( child.name( ) );
        PhysicalQuantity *quantity;

        if( name == PoPI_doubleChars ) {
            quantity = new PQ_double( child ); }
        else if( name == PoPI_integerChars ) {
            quantity = new PQ_integer( child ); }
        else if( name == PoPI_fractionChars ) {
            quantity = new PQ_fraction( child ); }
        else if( name == PoPI_stringChars ) {
            quantity = new PQ_string( child ); }
        else if( name == PoPI_shellChars ) {
            quantity = new PQ_shell( child ); }
        else {
            continue;
        }
        push_back( quantity );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

PQ_suite::~PQ_suite( ) {

    std::string::size_type i1, __size = size( );

    for( i1 = 0; i1 < __size; ++i1 ) delete (*this)[i1];
}

/* *********************************************************************************************************//**
 * Adds the contents of *this* to *a_XMLList* where each item in *a_XMLList* is one line (without linefeeds) to output as an XML representation of *this*.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void PQ_suite::toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    std::string indent2 = a_indent1 + "  ";

    if( size( ) == 0 ) return;
    std::string header = a_indent1 + "<" + m_label + ">";
    a_XMLList.push_back( std::move( header ) );
    for( std::vector<PhysicalQuantity *>::const_iterator iter = begin( ); iter != end( ); ++iter )
        (*iter)->toXMLList( a_XMLList, indent2 );
    appendXMLEnd( a_XMLList, m_label );
}

}
