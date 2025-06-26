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

#define PoPI_A_Chars "A"
#define PoPI_isotopeChars "isotope"

/*! \class Isotope
 * This class represents **PoPs** isotope instance.
 */

/* *********************************************************************************************************//**
 * Constructor that parses an **HAPI** instance to create a **PoPs** isotope node.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_DB              [in]    The **PoPI::Database:: instance to add the constructed **Isotope** to.
 * @param a_chemicalElement [in]    The parent chemical element suite that will contain *this*.
 ***********************************************************************************************************/

Isotope::Isotope( HAPI::Node const &a_node, Database *a_DB, ChemicalElement *a_chemicalElement ) :
        SymbolBase( a_node, Particle_class::isotope ),
        m_chemicalElement( a_chemicalElement ),
        m_Z( a_chemicalElement->Z( ) ),
        m_A( a_node.attribute( PoPI_A_Chars ).as_int( ) ),
        m_nuclides( PoPI_nuclidesChars ) {

    m_nuclides.appendFromParentNode( a_node.child( PoPI_nuclidesChars ), a_DB, this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Isotope::~Isotope( ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

void Isotope::calculateNuclideGammaBranchStateInfos( PoPI::Database const &a_pops, NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos ) const {

    NuclideGammaBranchStateInfo *info = nullptr, *priorInfo = nullptr;
    std::string energyUnit = "MeV";                     // FIXME, MeV needs to be a variable.
    double width = 0.0;

    for( std::size_t i1 = 0; i1 <  m_nuclides.size( ); ++i1 ) {
        std::size_t initialSize = a_nuclideGammaBranchStateInfos.size( );
        Nuclide const &nuclide = m_nuclides[i1];
        double levelEnergy = nuclide.levelEnergy( energyUnit );

        nuclide.calculateNuclideGammaBranchStateInfos( a_pops, a_nuclideGammaBranchStateInfos );
        if( a_nuclideGammaBranchStateInfos.size( ) > initialSize ) {
            info = a_nuclideGammaBranchStateInfos[initialSize];
            if( priorInfo != nullptr ) {
                if( priorInfo->kind( ) == PoPI_continuumChars ) {
                    for( std::size_t i2 = i1 + 1; i2 <  m_nuclides.size( ); ++i2 ) {
                        Nuclide const &nuclide2 = m_nuclides[i2];
                        if( nuclide2.levelEnergy( energyUnit ) > levelEnergy ) {
                            width = nuclide2.levelEnergy( energyUnit ) - levelEnergy;
                            break;
                        }
                    }
                    priorInfo->setNuclearLevelEnergyWidth( width );
                }
            }
        }

        priorInfo = info;
    }
}

/* *********************************************************************************************************//**
 * Adds the contents of *this* to *a_XMLList* where each item in *a_XMLList* is one line (without linefeeds) to output as an XML representation of *this*.
 *
 * @param a_XMLList                     [in]    The list to add an XML output representation of *this* to.
 * @param a_indent1                     [in]    The amount of indentation to added to each line added to *a_XMLList*.
 ***********************************************************************************************************/

void Isotope::toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    std::string::size_type size = m_nuclides.size( );
    std::string AStr = LUPI::Misc::argumentsToString( "%d", m_A );

    std::string header = a_indent1 + "<isotope symbol=\"" + symbol( ) + "\" A=\"" + AStr + "\">";
    a_XMLList.push_back( header );

    std::string indent2 = a_indent1 + "  ";
    std::string nuclideSuite = indent2 + "<" + PoPI_nuclidesChars + ">";
    a_XMLList.push_back( nuclideSuite );

    std::string indent3 = indent2 + "  ";
    for( std::string::size_type i1 = 0; i1 < size; ++i1 ) m_nuclides[i1].toXMLList( a_XMLList, indent3 );

    appendXMLEnd( a_XMLList, PoPI_nuclidesChars );
    appendXMLEnd( a_XMLList, PoPI_isotopeChars );
}

}
