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

#define PoPI_decayModesChars "decayModes"
#define PoPI_decayModeChars "decayMode"
#define PoPI_decayPathChars "decayPath"
#define PoPI_decayChars "decay"
#define PoPI_productsChars "products"
#define PoPI_photonEmissionProbabilitiesChars "photonEmissionProbabilities"

#define PoPI_typeChars "type"
#define PoPI_modeChars "mode"
#define PoPI_completeChars "complete"
#define PoPI_probabilityChars "probability"

/*
============================================================
======================== DecayData =========================
============================================================
*/
DecayData::DecayData( HAPI::Node const &a_node ) :
        m_decayModes( PoPI_decayModesChars ) {

    m_decayModes.appendFromParentNode2( a_node.child( PoPI_decayModesChars ), this );
}
/*
=========================================================
*/
DecayData::~DecayData( ) {

}
/*
=========================================================
*/
void DecayData::calculateNuclideGammaBranchStateInfo( PoPI::Database const &a_pops, NuclideGammaBranchStateInfo &a_nuclideGammaBranchStateInfo ) const {

    for( std::size_t i1 = 0; i1 <  m_decayModes.size( ); ++i1 ) {
        DecayMode const &decayMode = m_decayModes[i1];

        decayMode.calculateNuclideGammaBranchStateInfo( a_pops, a_nuclideGammaBranchStateInfo );
    }
}

/* *********************************************************************************************************//**
 * Adds the contents of *this* to *a_XMLList* where each item in *a_XMLList* is one line (without linefeeds) to output as an XML representation of *this*.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void DecayData::toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    std::string::size_type size = m_decayModes.size( );

    if( size == 0 ) return;

    std::string header = a_indent1 + "<" + PoPI_decayDataChars + ">";
    a_XMLList.push_back( std::move( header ) );

    if( size > 0 ) {
        std::string indent2 = a_indent1 + "  ";
        m_decayModes.toXMLList( a_XMLList, indent2 );
    }

    appendXMLEnd( a_XMLList, PoPI_decayDataChars );
}

/*
============================================================
======================== DecayMode =========================
============================================================
*/
DecayMode::DecayMode( HAPI::Node const &a_node, LUPI_maybeUnused DecayData const *a_decayData ) :
        m_label( a_node.attribute( PoPI_labelChars ).value( ) ),
        m_mode( a_node.attribute( PoPI_modeChars ).value( ) ),
        m_probability( a_node.child( PoPI_probabilityChars ) ),
        m_photonEmissionProbabilities( a_node.child( PoPI_photonEmissionProbabilitiesChars ) ),
        m_decayPath( PoPI_decayPathChars ) {

    m_decayPath.appendFromParentNode2( a_node.child( PoPI_decayPathChars ), this );
}
/*
============================================================
*/
DecayMode::~DecayMode( ) {

}

/*
============================================================
*/
void DecayMode::calculateNuclideGammaBranchStateInfo( PoPI::Database const &a_pops, NuclideGammaBranchStateInfo &a_nuclideGammaBranchStateInfo ) const {

    if( m_mode == PoPI_decayModeElectroMagnetic ) {
        double _probability = getPhysicalQuantityOfSuiteAsDouble( probability( ) );
        double _photonEmissionProbabilities = getPhysicalQuantityOfSuiteAsDouble( photonEmissionProbabilities( ), true, 1.0 );

        std::string residualState( "" );
        Decay const &decay = m_decayPath[0];
        Suite<Product, Decay> const &products = decay.products( );
        for( std::size_t i1 = 0; i1 < products.size( ); ++i1 ) {
            Product const &product = products[i1];

            if( product.pid( ) != IDs::photon ) residualState = product.pid( );
        }

        Particle const &initialState = a_pops.get<Particle>( a_nuclideGammaBranchStateInfo.state( ) );
        Particle const &finalState = a_pops.get<Particle>( residualState );
        double gammaEnergy = PoPI_AMU2MeV_c2 * ( initialState.massValue( "amu" ) - finalState.massValue( "amu" ) );

        NuclideGammaBranchInfo nuclideGammaBranchInfo( _probability, _photonEmissionProbabilities, gammaEnergy, residualState );
        a_nuclideGammaBranchStateInfo.add( nuclideGammaBranchInfo );
    }
}

/* *********************************************************************************************************//**
 * Adds the contents of *this* to *a_XMLList* where each item in *a_XMLList* is one line (without linefeeds) to output as an XML representation of *this*.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void DecayMode::toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    std::string header = a_indent1 + "<decayMode label=\"" + m_label + "\" mode=\"" + m_mode + "\">";
    a_XMLList.push_back( std::move( header ) );

    std::string indent2 = a_indent1 + "  ";
    m_probability.toXMLList( a_XMLList, indent2 );
    m_decayPath.toXMLList( a_XMLList, indent2 );

    appendXMLEnd( a_XMLList, PoPI_decayModeChars );
}

/*
============================================================
========================== Decay ===========================
============================================================
*/

Decay::Decay( HAPI::Node const &a_node, LUPI_maybeUnused DecayMode const *a_decayMode ) :
        m_index( a_node.attribute( PoPI_indexChars ).as_int( ) ),
        m_mode( a_node.attribute( PoPI_modeChars ).value( ) ),
        m_complete( a_node.attribute( PoPI_completeChars ).value( ) == "true" ),
        m_products( PoPI_productsChars ) {

    if( a_node.attribute( PoPI_typeChars ).value( ) != "" ) m_mode = a_node.attribute( PoPI_typeChars ).value( );

    m_products.appendFromParentNode2( a_node.child( PoPI_productsChars ), this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Decay::~Decay( ) {

}

/* *********************************************************************************************************//**
 * Adds the contents of *this* to *a_XMLList* where each item in *a_XMLList* is one line (without linefeeds) to output as an XML representation of *this*.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void Decay::toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    std::string indexString( std::to_string( m_index ) );

    std::string header = a_indent1 + "<decay index=\"" + indexString + "\"";
    if( m_mode != "" ) header += " mode=\"" + m_mode + "\"";
    if( m_complete ) header += " mode=\"true\"";
    header += ">";
    a_XMLList.push_back( std::move( header ) );

    std::string indent2 = a_indent1 + "  ";
    m_products.toXMLList( a_XMLList, indent2 );

    appendXMLEnd( a_XMLList, PoPI_decayChars );
}

/*
============================================================
========================= Product ==========================
============================================================
*/
Product::Product( HAPI::Node const &a_node, LUPI_maybeUnused Decay *a_DB ) :
        m_id( -1 ),
        m_pid( a_node.attribute( PoPI_pidChars ).value( ) ),
        m_label( a_node.attribute( PoPI_labelChars ).value( ) ) {

}

/*
============================================================
*/
Product::~Product( ) {

}

/* *********************************************************************************************************//**
 * Adds the contents of *this* to *a_XMLList* where each item in *a_XMLList* is one line (without linefeeds) to output as an XML representation of *this*.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void Product::toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    std::string header = a_indent1 + "<product label=\"" + m_label + "\" pid=\"" + m_pid + "\"/>";
    a_XMLList.push_back( std::move( header ) );
}

/*! \class GammaDecayData
 * This class stores, in a crude way, the GRIN, non-GNDS 2.0 compliant, nuclide gamma decay data.
 */

/* *********************************************************************************************************//**
 * @param a_node                        [in]    The HAPI node to parse.
 ***********************************************************************************************************/

GammaDecayData::GammaDecayData( HAPI::Node const &a_node ) :
        m_kind( a_node.attribute_as_string( "kind" ) ),
        m_rows( 0 ),
        m_columns( 0 ) {

    if( m_kind == "" ) m_kind = PoPI_discreteChars;

    if( !a_node.empty( ) ) {
        HAPI::Node table = a_node.child( "table" );
        m_rows = table.attribute_as_int( "rows" );
        m_columns = table.attribute_as_int( "columns" );

        HAPI::Node data = table.child( "data" );

        std::string text = LUPI::Misc::stripString( data.text( ).get( ) );
        auto cells = LUPI::Misc::splitString( text, ' ', true );

        m_ids.reserve( static_cast<std::size_t>( m_rows ) );
        m_probabilities.reserve( static_cast<std::size_t>( m_rows ) );
        m_photonEmissionProbabilities.reserve( static_cast<std::size_t>( m_rows ) );
        for( std::size_t cellIndex = 0; cellIndex < cells.size( ); cellIndex += 3 ) {
            m_ids.push_back( cells[cellIndex] );
            m_probabilities.push_back( std::stod( cells[cellIndex+1] ) );
            m_photonEmissionProbabilities.push_back( std::stod( cells[cellIndex+2] ) );
        }
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

GammaDecayData::~GammaDecayData( ) {

}

/*
============================================================
*/
void GammaDecayData::calculateNuclideGammaBranchStateInfo( PoPI::Database const &a_pops, NuclideGammaBranchStateInfo &a_nuclideGammaBranchStateInfo ) const {

    Particle const &initialState = a_pops.get<Particle>( a_nuclideGammaBranchStateInfo.state( ) );
    double initialStateMass = initialState.massValue( "amu" );

    for( std::size_t index = 0; index < static_cast<std::size_t>( m_rows ); ++index ) {
        std::string residualState( m_ids[index] );
        double _probability = m_probabilities[index];
        double _photonEmissionProbabilities = m_photonEmissionProbabilities[index];

        Particle const &finalState = a_pops.get<Particle>( residualState );
        double gammaEnergy = PoPI_AMU2MeV_c2 * ( initialStateMass - finalState.massValue( "amu" ) );

        NuclideGammaBranchInfo nuclideGammaBranchInfo( _probability, _photonEmissionProbabilities, gammaEnergy, residualState );
        a_nuclideGammaBranchStateInfo.add( nuclideGammaBranchInfo );
    }
}

}
