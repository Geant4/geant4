/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "GIDI.hpp"

namespace GIDI {

namespace ACE_URR {

/*! \class ProbabilityTable
 * Class for the LLNL defined **probabilityTable** node which is not a part of the **GNDS* 2.0 specifications. 
 * This node currently can be found in the **applicationData** node of a **reactionSuite** file with the label 
 * "*LLNL::URR_probability_tables*".
 */

/* *********************************************************************************************************//**
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

ProbabilityTable::ProbabilityTable( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Form( a_node, a_setupInfo, FormType::ACE_URR_probabilityTable, a_parent ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        m_forms.push_back( new IncidentEnergy( a_construction, child, a_setupInfo ) );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ProbabilityTable::~ProbabilityTable( ) {

    for( auto iter = m_forms.begin( ); iter != m_forms.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void ProbabilityTable::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes;

    attributes = a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    for( auto iter = m_forms.begin( ); iter != m_forms.end( ); ++iter ) (*iter)->toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class IncidentEnergy
 * Class for the LLNL define **incidentEnergy** node which is not a **GNDS* 2.0 specifications. This node currently can be
 * found in the **applicationData** node of a **reactionSuite** file with the label "*LLNL::ACE_URR_probability_tables*".
 */

/* *********************************************************************************************************//**
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

IncidentEnergy::IncidentEnergy( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::ACE_URR_probabilityTable ),
        m_value( a_node.attribute( GIDI_valueChars ).as_double( ) ),
        m_unit( a_node.attribute_as_string( GIDI_unitChars ) ),
        m_table( a_construction, a_node.first_child( ), a_setupInfo ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

IncidentEnergy::~IncidentEnergy( ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void IncidentEnergy::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes;

    attributes = a_writeInfo.addAttribute( GIDI_valueChars, LUPI::Misc::doubleToShortestString( m_value ) );
    attributes += a_writeInfo.addAttribute( GIDI_unitChars, m_unit );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    m_table.toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeEnder( moniker( ) );
}

}               // End namespace ACE_URR.

}               // End namespace GIDI.
