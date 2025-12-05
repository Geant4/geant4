/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <algorithm>
#include <cmath>

#include <GIDI.hpp>
#include <HAPI.hpp>

namespace GIDI {

namespace TargetInfo {

static GUPI::Entry *parseChemicalElement( LUPI_maybeUnused GUPI::Suite *a_parent, HAPI::Node const &a_node );
static GUPI::Entry *parseNuclide( LUPI_maybeUnused GUPI::Suite *a_parent, HAPI::Node const &a_node );

/*! \class Nuclide
 * The class that stores a chemicalElement node.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Nuclide::Nuclide( HAPI::Node const &a_node ) :
        GUPI::Entry( a_node, GIDI_pidChars ),
        m_atomFraction( a_node.attribute( GIDI_atomFractionChars ).as_double( ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Nuclide::~Nuclide( ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Nuclide::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;

    attributes  = a_writeInfo.addAttribute( GIDI_pidChars, pid( ) );
    attributes += a_writeInfo.addAttribute( GIDI_atomFractionChars, LUPI::Misc::doubleToShortestString( atomFraction( ), 15, -2 ) );
    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}

/*! \class ChemicalElement
 * The class that stores a chemicalElement node.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ChemicalElement::ChemicalElement( HAPI::Node const &a_node ) :
        GUPI::Entry( a_node, PoPI_symbolChars ),
        m_nuclides( a_node.child( PoPI_nuclidesChars ), GIDI_pidChars, parseNuclide ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ChemicalElement::~ChemicalElement( ) {

}

/* *********************************************************************************************************//**
 * Returns the **Nuclide** with pid *a_pid* if it exists; otherwise, **nullptr** is returned.
 *
 * @param a_symbol      [in]    The symbol for the chemical element whose nuclide abundance data are being requested.
 ***********************************************************************************************************/

Nuclide const *ChemicalElement::operator[]( std::string const &a_pid ) const {

    for( auto nuclideIter = m_nuclides.begin( ); nuclideIter != m_nuclides.end( ); ++nuclideIter ) {
        Nuclide const *nuclide = static_cast<Nuclide *>( *nuclideIter );
        if( nuclide->pid( ) == a_pid ) return( nuclide );
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void ChemicalElement::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes;

    attributes  = a_writeInfo.addAttribute( PoPI_symbolChars, symbol( ) );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    m_nuclides.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class IsotopicAbundances
 * The class that stores an isotopicAbundances node.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

IsotopicAbundances::IsotopicAbundances( ) :
        GUPI::Ancestry( GIDI_isotopicAbundancesChars ),
        m_chemicalElements( PoPI_chemicalElementsChars, PoPI_symbolChars ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

IsotopicAbundances::~IsotopicAbundances( ) {

}

void IsotopicAbundances::initialize( HAPI::Node const &a_node ) {

    m_chemicalElements.parse( a_node.child( PoPI_chemicalElementsChars ), parseChemicalElement );
}

/* *********************************************************************************************************//**
 * Returns the **ChemicalElement** with symbol *a_symbol* if it exists; otherwise, **nullptr** is returned.
 *
 * @param a_symbol      [in]    The symbol for the chemical element whose isotopic abundance data are being requested.
 ***********************************************************************************************************/

ChemicalElement const *IsotopicAbundances::operator[]( std::string const &a_symbol ) const {

    for( auto iter = m_chemicalElements.begin( ); iter != m_chemicalElements.end( ); ++iter ) {
        ChemicalElement const *chemicalElement = static_cast<ChemicalElement const *>( *iter );
        if( chemicalElement->symbol( ) == a_symbol ) return( chemicalElement );
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void IsotopicAbundances::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ) );
    m_chemicalElements.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class TargetInfo
 * The class that stores the **targetInfo** node.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

TargetInfo::TargetInfo( ) :
        GUPI::Ancestry( GIDI_targetInfoChars ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

TargetInfo::~TargetInfo( ) {

}

/* *********************************************************************************************************//**
 * For internal use only.
 *
 * @param a_node                [in]    The **HAPI::Node** node whose text is to be converted into a list of doubles.
 ***********************************************************************************************************/

void TargetInfo::parseEvaluatedTargetInfo( HAPI::Node const &a_node ) {

    m_isotopicAbundances.initialize( a_node.child( GIDI_isotopicAbundancesChars ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void TargetInfo::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ) );
    m_isotopicAbundances.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

static GUPI::Entry *parseChemicalElement( LUPI_maybeUnused GUPI::Suite *a_parent, HAPI::Node const &a_node ) {

    return new ChemicalElement( a_node );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

static GUPI::Entry *parseNuclide( LUPI_maybeUnused GUPI::Suite *a_parent, HAPI::Node const &a_node ) {

    return new Nuclide( a_node );
}

}               // End of namespace TargetInfo.

}               // End of namespace GIDI.
