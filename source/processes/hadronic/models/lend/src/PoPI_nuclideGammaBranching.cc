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

/*! \class NuclideGammaBranchInfo
 * Class storing information about the gamma (i.e., photon) decay of an excited nuclide state to a lower state.
 */

/* *********************************************************************************************************//**
 * @param a_probability                 [in]    The probability that the level decays to state *a_residualState*.
 * @param a_photonEmissionProbability   [in]    The conditional probability the the decay emitted a photon.
 * @param a_gammaEnergy                 [in]    The energy of the emitted photon.
 * @param a_residualState               [in]    The state the residual is left in after photon decay.
 ***********************************************************************************************************/

NuclideGammaBranchInfo::NuclideGammaBranchInfo( double a_probability, double a_photonEmissionProbability, double a_gammaEnergy, 
                std::string const &a_residualState ) :
        m_probability( a_probability ),
        m_photonEmissionProbability( a_photonEmissionProbability ),
        m_gammaEnergy( a_gammaEnergy ),
        m_residualState( a_residualState ) {

}

/* *********************************************************************************************************//**
 * Copy constructor.

 * @param a_nuclideGammaBranchInfo      [in]    The *NuclideGammaBranchInfo* instance to copy.
 ***********************************************************************************************************/

NuclideGammaBranchInfo::NuclideGammaBranchInfo( NuclideGammaBranchInfo const &a_nuclideGammaBranchInfo ) :
        m_probability( a_nuclideGammaBranchInfo.probability( ) ),
        m_photonEmissionProbability( a_nuclideGammaBranchInfo.photonEmissionProbability( ) ),
        m_gammaEnergy( a_nuclideGammaBranchInfo.gammaEnergy( ) ),
        m_residualState( a_nuclideGammaBranchInfo.residualState( ) ) {

}

/*
============================================================
================= NuclideGammaBranchStateInfo ================
============================================================
*/
NuclideGammaBranchStateInfo::NuclideGammaBranchStateInfo( std::string const &a_state, int a_intid, std::string const &a_kind, 
                double a_nuclearLevelEnergy ) :
        m_state( a_state ),
        m_intid( a_intid ),
        m_kind( a_kind ),
        m_nuclearLevelEnergy( a_nuclearLevelEnergy ),
        m_nuclearLevelEnergyWidth( 0.0 ),
        m_derivedCalculated( false ),
        m_multiplicity( 0.0 ),
        m_averageGammaEnergy( 0.0 ) {

}
/*
=========================================================
*/
void NuclideGammaBranchStateInfo::add( NuclideGammaBranchInfo const &a_nuclideGammaBranchInfo ) {

    m_branches.push_back( a_nuclideGammaBranchInfo );
}
/*
=========================================================
*/
void NuclideGammaBranchStateInfo::calculateDerivedData( NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos ) {

    if( m_derivedCalculated ) return;

    for( std::size_t i1 = 0; i1 < m_branches.size( ); ++i1 ) {
        NuclideGammaBranchInfo &nuclideGammaBranchInfo = m_branches[i1];

        std::string const &residualState = nuclideGammaBranchInfo.residualState( );
        NuclideGammaBranchStateInfo *nuclideGammaBranchStateInfo = a_nuclideGammaBranchStateInfos.find( residualState );

        double chainedMultiplicity = 0.0;
        double chainedAverageGammaEnergy = 0.0;
        if( nuclideGammaBranchStateInfo != nullptr ) {
            nuclideGammaBranchStateInfo->calculateDerivedData( a_nuclideGammaBranchStateInfos );
            chainedMultiplicity = nuclideGammaBranchStateInfo->multiplicity( );
            chainedAverageGammaEnergy = nuclideGammaBranchStateInfo->averageGammaEnergy( );
        }

        m_multiplicity += nuclideGammaBranchInfo.probability( ) * ( nuclideGammaBranchInfo.photonEmissionProbability( ) + chainedMultiplicity );
        m_averageGammaEnergy += nuclideGammaBranchInfo.probability( ) * 
                ( nuclideGammaBranchInfo.photonEmissionProbability( ) * nuclideGammaBranchInfo.gammaEnergy( ) + chainedAverageGammaEnergy );
    }

    m_derivedCalculated = true;
}

/*
============================================================
================ NuclideGammaBranchStateInfos ================
============================================================
*/
NuclideGammaBranchStateInfos::NuclideGammaBranchStateInfos( ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

NuclideGammaBranchStateInfos::~NuclideGammaBranchStateInfos( ) {

    for( std::size_t i1 = 0; i1 < m_nuclideGammaBranchStateInfos.size( ); ++i1 ) delete m_nuclideGammaBranchStateInfos[i1];
}
/*
=========================================================
*/
void NuclideGammaBranchStateInfos::add( NuclideGammaBranchStateInfo *a_nuclideGammaBranchStateInfo ) {

    m_nuclideGammaBranchStateInfos.push_back( a_nuclideGammaBranchStateInfo );
}

/* *********************************************************************************************************//**
 * This method returns a pointer to the NuclideGammaBranchStateInfo instance for *a_state* or nullptr if not match is found.
 *
 * @param a_state           [in]    The PoPs id for the requested state (i.e., nuclide).
 *
 * @return                          A pointer to the requested NuclideGammaBranchStateInfo instance or nullptr if not match is found.
 ***********************************************************************************************************/

NuclideGammaBranchStateInfo *NuclideGammaBranchStateInfos::find( std::string const &a_state ) {

    for( std::size_t i1 = 0; i1 < m_nuclideGammaBranchStateInfos.size( ); ++i1 ) {
        NuclideGammaBranchStateInfo *nuclideGammaBranchStateInfo = m_nuclideGammaBranchStateInfos[i1];

        if( nuclideGammaBranchStateInfo->state( ) == a_state ) return( nuclideGammaBranchStateInfo );
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * This method returns a const pointer to the NuclideGammaBranchStateInfo instance for *a_state* or nullptr if not match is found.
 *
 * @param a_state           [in]    The PoPs id for the requested state (i.e., nuclide).
 *
 * @return                          A const pointer to the requested NuclideGammaBranchStateInfo instance or nullptr if not match is found.
 ***********************************************************************************************************/

NuclideGammaBranchStateInfo const *NuclideGammaBranchStateInfos::find( std::string const &a_state ) const {

    for( std::size_t i1 = 0; i1 < m_nuclideGammaBranchStateInfos.size( ); ++i1 ) {
        NuclideGammaBranchStateInfo *nuclideGammaBranchStateInfo = m_nuclideGammaBranchStateInfos[i1];

        if( nuclideGammaBranchStateInfo->state( ) == a_state ) return( nuclideGammaBranchStateInfo );
    }

    return( nullptr );
}

}
