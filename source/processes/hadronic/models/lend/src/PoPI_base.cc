/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "PoPI.hpp"

#define PoPI_idChars "id"
#define PoPI_symbolChars "symbol"

namespace PoPI {

/*! \class Base
 * This class is the base class for all **Particle** and **SymbolBase** instances.
 */

/* *********************************************************************************************************//**
 * @param a_id              [in]    The **PoPs** id for *this*.
 * @param a_class           [in]    The **PoPI** class for *this*.
 ***********************************************************************************************************/

Base::Base( std::string const &a_id, Particle_class a_class ) :
        m_id( a_id ),
        m_class( a_class ),
        m_index( -1 ),
        m_intid( -1 ) {

}

/* *********************************************************************************************************//**
 * Constructor that parses an **HAPI** instance.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_label           [in]    This is either *id* or *symbol*. That is, it is the name of the attribute in *a_node* whose value represents *this* *m_id* member.
 * @param a_class           [in]    The **PoPI** class for *this*.
 ***********************************************************************************************************/

Base::Base( HAPI::Node const &a_node, std::string const &a_label, Particle_class a_class ) :
        m_id( a_node.attribute( a_label.c_str( ) ).value( ) ),
        m_class( a_class ),
        m_index( -1 ),
        m_intid( -1 ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Base::~Base( ) {

}

/*! \class IDBase
 * This class is the base class for all **Particle** instances.
 */

/* *********************************************************************************************************//**
 * @param a_id              [in]    The **PoPs** id for *this*.
 * @param a_class           [in]    The **PoPI** class for *this*.
 ***********************************************************************************************************/

IDBase::IDBase( std::string const &a_id, Particle_class a_class ) :
        Base( a_id, a_class ) {

}

/* *********************************************************************************************************//**
 * Constructor that parses an **HAPI** instance.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_class           [in]    The **PoPI** class for *this*.
 ***********************************************************************************************************/

IDBase::IDBase( HAPI::Node const &a_node, Particle_class a_class ) :
        Base( a_node, PoPI_idChars, a_class ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

IDBase::~IDBase( ) {

}

/* *********************************************************************************************************//**
 * This method adds *this* to the *m_list* member of *a_DB*.
 *
 * @param a_DB              [in]    The **PoPI::Database** instance to add *this* to.
 *
 * @return                          The index assigned to *this* by *a_DB*.
 ***********************************************************************************************************/

int IDBase::addToDatabase( Database *a_DB ) {

    a_DB->add( this );
    return( index( ) );
}

/*! \class SymbolBase
 * This class is the base class for all **SymbolBase** instances.
 */

/* *********************************************************************************************************//**
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_class           [in]    The **PoPI** class for *this*.
 ***********************************************************************************************************/

SymbolBase::SymbolBase( HAPI::Node const &a_node, Particle_class a_class ) :
        Base( a_node, PoPI_symbolChars, a_class ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

SymbolBase::~SymbolBase( ) {

}

/* *********************************************************************************************************//**
 * This method adds *this* to the *m_symbolList* member of *a_DB*.
 *
 * @param a_DB              [in]    The **PoPI::Database** instance to add *this* to.
 *
 * @return                          The index assigned to *this* by *a_DB*.
 ***********************************************************************************************************/

int SymbolBase::addToSymbols( Database *a_DB ) {

    a_DB->addSymbol( this );
    return( index( ) );
}

}
