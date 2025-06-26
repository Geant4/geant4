/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <string.h>
#include <iomanip>

#include "MCGIDI.hpp"

namespace MCGIDI {

/*
============================================================
================== NuclideGammaBranchInfo ==================
============================================================
*/

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE NuclideGammaBranchInfo::NuclideGammaBranchInfo( ) :
        m_probability( 0.0 ),
        m_photonEmissionProbability( 0.0 ),
        m_gammaEnergy( 0.0 ),
        m_residualStateIndex( -1 ),
        m_residualStateKindIsContinuum( false ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

NuclideGammaBranchInfo::NuclideGammaBranchInfo( PoPI::NuclideGammaBranchInfo const &a_nuclideGammaBranchInfo, 
                std::map<std::string, int> &a_stateNamesToIndices, bool a_makePhotonEmissionProbabilitiesOne ) :
        m_probability( a_nuclideGammaBranchInfo.probability( ) ),
        m_photonEmissionProbability( a_nuclideGammaBranchInfo.photonEmissionProbability( ) ),
        m_gammaEnergy( a_nuclideGammaBranchInfo.gammaEnergy( ) ),
        m_residualStateIndex( -1 ),
        m_residualStateKindIsContinuum( false ) {

        if( a_makePhotonEmissionProbabilitiesOne ) m_photonEmissionProbability = 1.0;
        std::map<std::string,int>::iterator iter = a_stateNamesToIndices.find( a_nuclideGammaBranchInfo.residualState( ) );
        if( iter != a_stateNamesToIndices.end( ) ) m_residualStateIndex = iter->second;
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void NuclideGammaBranchInfo::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_DOUBLE( m_probability, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_photonEmissionProbability, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_gammaEnergy, a_buffer, a_mode );
    DATA_MEMBER_INT( m_residualStateIndex, a_buffer, a_mode );
    DATA_MEMBER_CAST( m_residualStateKindIsContinuum, a_buffer, a_mode, bool );
}

/* *********************************************************************************************************//**
 * Print to *std::cout* the content of *this*. This is mainly meant for debugging.
 *
 * @param a_protareSingle       [in]    The ProtareSingle instance *this* resides in.
 * @param a_indent              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_iFormat             [in]    C printf format specifier for any interger that is printed (e.g., "%3d").
 * @param a_energyFormat        [in]    C printf format specifier for any interger that is printed (e.g., "%20.12e").
 * @param a_dFormat             [in]    C printf format specifier for any interger that is printed (e.g., "%14.7e").
 ***********************************************************************************************************/
LUPI_HOST void NuclideGammaBranchInfo::print( LUPI_maybeUnused ProtareSingle const *a_protareSingle, std::string const &a_indent, std::string const &a_iFormat,
                LUPI_maybeUnused std::string const &a_energyFormat, std::string const &a_dFormat ) const {

    std::cout << a_indent << std::left << std::setw( 17 ) << LUPI::Misc::argumentsToString( a_dFormat.c_str( ), m_probability )
            << LUPI::Misc::argumentsToString( a_dFormat.c_str( ), m_photonEmissionProbability )
            << LUPI::Misc::argumentsToString( a_dFormat.c_str( ), m_gammaEnergy )
            << LUPI::Misc::argumentsToString( a_iFormat.c_str( ), m_residualStateIndex ) << std::endl;
}

/*
============================================================
============== NuclideGammaBranchStateInfo =================
============================================================
*/

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE NuclideGammaBranchStateInfo::NuclideGammaBranchStateInfo( ) :
        m_intid( -1 ),
        m_nuclearLevelEnergy( 0.0 ),
        m_nuclearLevelEnergyWidth( 0.0 ),
        m_multiplicity( 0.0 ),
        m_averageGammaEnergy( 0.0 ) {

        m_state[0] = 0;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

NuclideGammaBranchStateInfo::NuclideGammaBranchStateInfo( PoPI::NuclideGammaBranchStateInfo const &a_nuclideGammaBranchingInfo,
                std::vector<NuclideGammaBranchInfo *> &a_nuclideGammaBranchInfos, 
                std::map<std::string, int> &a_stateNamesToIndices, bool a_makePhotonEmissionProbabilitiesOne,
                bool a_zeroNuclearLevelEnergyWidth ) :
        m_intid( a_nuclideGammaBranchingInfo.intid( ) ),
        m_nuclearLevelEnergy( a_nuclideGammaBranchingInfo.nuclearLevelEnergy( ) ),
        m_nuclearLevelEnergyWidth( a_nuclideGammaBranchingInfo.nuclearLevelEnergyWidth( ) ),
        m_multiplicity( a_nuclideGammaBranchingInfo.multiplicity( ) ),
        m_averageGammaEnergy( a_nuclideGammaBranchingInfo.averageGammaEnergy( ) ) {

    if( a_zeroNuclearLevelEnergyWidth ) m_nuclearLevelEnergyWidth = 0.0;

    strncpy( m_state, a_nuclideGammaBranchingInfo.state( ).c_str( ), sizeof( m_state ) );
    m_state[sizeof( m_state )-1] = 0;

    std::vector<PoPI::NuclideGammaBranchInfo> const &branches = a_nuclideGammaBranchingInfo.branches( );
    m_branchIndices.reserve( branches.size( ) );

    for( std::size_t i1 = 0; i1 < branches.size( ); ++i1 ) {
        m_branchIndices.push_back( a_nuclideGammaBranchInfos.size( ) );
        a_nuclideGammaBranchInfos.push_back( new NuclideGammaBranchInfo( branches[i1], a_stateNamesToIndices, a_makePhotonEmissionProbabilitiesOne ) );
    }
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void NuclideGammaBranchStateInfo::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_CHAR_ARRAY( m_state, a_buffer, a_mode );
    DATA_MEMBER_INT( m_intid, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_nuclearLevelEnergy, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_nuclearLevelEnergyWidth, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_multiplicity, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_averageGammaEnergy, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_branchIndices, a_buffer, a_mode );
}

/* *********************************************************************************************************//**
 * Print to *std::cout* the content of *this*. This is mainly meant for debugging.
 *
 * @param a_protareSingle       [in]    The ProtareSingle instance *this* resides in.
 * @param a_indent              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_iFormat             [in]    C printf format specifier for any interger that is printed (e.g., "%3d").
 * @param a_energyFormat        [in]    C printf format specifier for any interger that is printed (e.g., "%20.12e").
 * @param a_dFormat             [in]    C printf format specifier for any interger that is printed (e.g., "%14.7e").
 ***********************************************************************************************************/
LUPI_HOST void NuclideGammaBranchStateInfo::print( LUPI_maybeUnused ProtareSingle const *a_protareSingle, std::string const &a_indent, std::string const &a_iFormat,
                LUPI_maybeUnused std::string const &a_energyFormat, std::string const &a_dFormat ) const {

    std::cout << a_indent << std::left << std::setw( 17 ) << m_state << LUPI::Misc::argumentsToString( a_dFormat.c_str( ), m_multiplicity ) <<
            LUPI::Misc::argumentsToString( a_dFormat.c_str( ), m_averageGammaEnergy );
    for( auto branchIter = m_branchIndices.begin( ); branchIter != m_branchIndices.end( ); ++branchIter ) {
        std::cout << LUPI::Misc::argumentsToString( a_iFormat.c_str( ), (*branchIter) );
    }
    std::cout << std::endl;
}

}
