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

/*! \class Entry
 * This is a base class inherit by other classes that are enties in a **Suite**.
 */

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The **GNDS** node's name (i.e., moniker).
 * @param a_attribute           [in]    Currently not used.
 ***********************************************************************************************************/

Entry::Entry( std::string const &a_moniker, std::string const &a_keyName, std::string const &a_keyValue ) :
        Ancestry( a_moniker ),
        m_keyName( a_keyName ),
        m_keyValue( a_keyValue ) {

}

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The **GNDS** node's name (i.e., moniker).
 * @param a_attribute           [in]    Currently not used.
 ***********************************************************************************************************/

Entry::Entry( HAPI::Node const &a_node, std::string const &a_keyName ) :
        Ancestry( a_node.name( ) ),
        m_keyName( a_keyName ),
        m_keyValue( a_node.attribute_as_string( a_keyName.c_str( ) ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Entry::~Entry( ) {

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST void Entry::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Ancestry::serialize( a_buffer, a_mode );
    DATA_MEMBER_STD_STRING( m_keyName, a_buffer, a_mode );
    DATA_MEMBER_STD_STRING( m_keyValue, a_buffer, a_mode );
}

}
