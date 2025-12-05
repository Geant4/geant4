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

/*! \class PhysicalQuantity
 * Class to store a physical quantity. A physical quantity is a value (e.g., 13.2) with a unit (e.g., 'cm'). The physical quantity
 * can be unitless (i.e., the unit can be an empty string). Examples a physical quantities are '13.2 cm', '0.132 m', '4.5 kg'.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node            [in]     The **HAPI::Node** to be parsed and used to construct the PhysicalQuantity.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

PhysicalQuantity::PhysicalQuantity( HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::physicalQuantity ),
        m_value( a_node.attribute_as_double( GIDI_valueChars ) ),
        m_unit( a_node.attribute_as_string( GIDI_unitChars ) ) {

}

/* *********************************************************************************************************//**
 *
 * @param a_value           [in]     The physical quantity's value.
 * @param a_unit            [in]     The physical quantity's unit.
 ***********************************************************************************************************/

PhysicalQuantity::PhysicalQuantity( double a_value, std::string const &a_unit ) :
        Form( FormType::physicalQuantity ),
        m_value( a_value ),
        m_unit( a_unit ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

PhysicalQuantity::~PhysicalQuantity( ) {

}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the members of *this* to those of *a_rhs* except for those
 * not set by base classes.
 *
 * @param a_rhs                     [in]    Instance whose member are used to set the members of *this*.
 ***********************************************************************************************************/

PhysicalQuantity &PhysicalQuantity::operator=( PhysicalQuantity const &a_rhs ) {

    if( this != &a_rhs ) {
        Form::operator=( a_rhs );

        m_value = a_rhs.value( );
        m_unit = a_rhs.unit( );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void PhysicalQuantity::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes = a_writeInfo.addAttribute( GIDI_valueChars, LUPI::Misc::doubleToShortestString( value( ) ) ) + a_writeInfo.addAttribute( GIDI_unitChars, unit( ) );

    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}

/* *********************************************************************************************************//**
 * Writes the information of *a_physicalQuantity* to *a_os*.
 *
 * @param       a_os                [out]       The stream to write to.
 * @param       a_physicalQuantity  [in]        The PhysicalQuantity whose information is written.
 ***********************************************************************************************************/

std::ostream &operator<<( std::ostream &a_os, PhysicalQuantity const &a_physicalQuantity ) {

    a_os << a_physicalQuantity.value( ) << " " << a_physicalQuantity.unit( );

    return( a_os );
}

}
