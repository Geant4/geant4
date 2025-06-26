/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <MCGIDI.hpp>

namespace MCGIDI {

/*! \class GRIN_levelsAndProbabilities
 * This class stores a Vector of summed probabilities and a Vector of their associated nuclide level as needed by
 * inelastic and cpature GRIN continuum reaction data.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_levelsAndProbabilities::GRIN_levelsAndProbabilities( ) {

}

/* *********************************************************************************************************//**
 * @param a_setupInfo           [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param a_table               [in]    The table with a column containing nucide ids and a column with their probabilities.
 ***********************************************************************************************************/

LUPI_HOST GRIN_levelsAndProbabilities::GRIN_levelsAndProbabilities( SetupInfo &a_setupInfo, PoPI::Database const &a_pops, 
                GIDI::Table::Table const &a_table, bool a_normalize ) {

    m_summedProbabilities.reserve( a_table.rows( ) );
    m_levels.reserve( a_table.rows( ) );
    m_isModelledLevel.reserve( a_table.rows( ) );

    double sum = 0.0;
    auto cells = LUPI::Misc::splitString( LUPI::Misc::stripString( a_table.data( ).body( ) ), ' ', true );
    for( std::size_t index = 0; index < cells.size( ); ++index ) {
        m_levels.push_back( a_setupInfo.m_stateNamesToIndices[cells[index]] );
        PoPI::Nuclide const &nuclide = a_pops.get<PoPI::Nuclide const>( cells[index] );

        ++index;
        sum += std::stod( cells[index] );
        m_summedProbabilities.push_back( sum );

        m_isModelledLevel.push_back( nuclide.kind( ) == PoPI_continuumChars );
    }

    if( a_normalize ) {
        for( auto iter = m_summedProbabilities.begin( ); iter != m_summedProbabilities.end( ); ++iter ) (*iter) /= sum;
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_levelsAndProbabilities::~GRIN_levelsAndProbabilities( ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST void GRIN_levelsAndProbabilities::set( std::vector<int> const &a_levels, std::vector<double> const &a_probabilities ) {

    m_levels.reserve( a_levels.size( ) );
    m_summedProbabilities.reserve( a_levels.size( ) );
    m_isModelledLevel.reserve( a_levels.size( ) );

    double sum = 0;
    for( std::size_t index = 0; index < a_levels.size( ); ++index ) {
        m_levels.push_back( a_levels[index] );
        sum += a_probabilities[index];
        m_summedProbabilities.push_back( sum );
        m_isModelledLevel.push_back( true );
    }

    for( auto iter = m_summedProbabilities.begin( ); iter != m_summedProbabilities.end( ); ++iter ) {
        *iter /= sum;
    }
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void GRIN_levelsAndProbabilities::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_VECTOR_INT( m_levels, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_summedProbabilities, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_BOOL( m_isModelledLevel, a_buffer, a_mode );
}

/*! \class GRIN_inelasticForEnergy
 * This class represents GRIN inelastic continuum reaction data which has simulated levels.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_inelasticForEnergy::GRIN_inelasticForEnergy( ) {

}

/* *********************************************************************************************************//**
 * @param a_setupInfo                   [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param GRIN_continuumGammas          [in]    GIDI instance containing the GRIN capture data.
 ***********************************************************************************************************/

LUPI_HOST GRIN_inelasticForEnergy::GRIN_inelasticForEnergy( SetupInfo &a_setupInfo, double a_projectileMass, double a_targetMass, 
                PoPI::Database const &a_pops, GIDI::GRIN::InelasticIncidentEnergy const *inelasticIncidentEnergy ) :
        m_levelsAndProbabilities( a_setupInfo, a_pops, inelasticIncidentEnergy->table( ), true ) {

    std::vector<int> indices;
    std::vector<double> thresholds;
    int index = 0;
    double priorThreshold = -1;
    for( auto iter = m_levelsAndProbabilities.m_levels.begin( ); iter != m_levelsAndProbabilities.m_levels.end( ); ++iter, ++index ) {
        NuclideGammaBranchStateInfo const *nuclideGammaBranchStateInfo = a_setupInfo.m_protare.nuclideGammaBranchStateInfos( )[*iter];
        double levelEnergy = nuclideGammaBranchStateInfo->nuclearLevelEnergy( );

        double threshold = ( a_projectileMass + a_targetMass + levelEnergy / 2 ) * levelEnergy / a_targetMass;
        if( threshold > priorThreshold ) {
            if( thresholds.size( ) > 0 )
                indices.push_back( index - 1 );
            thresholds.push_back( threshold );
            priorThreshold = threshold;
        }
    }
    indices.push_back( m_levelsAndProbabilities.m_levels.size( ) - 1 );

    m_indices = indices;
    m_thresholds = thresholds;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_inelasticForEnergy::~GRIN_inelasticForEnergy( ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE int GRIN_inelasticForEnergy::sampleLevelIndex( double a_projectileEnergy, double a_random ) const {

    std::size_t index = 0;

    for( auto iter = m_thresholds.begin( ); iter != m_thresholds.end( ); ++iter, ++index ) {
        if( *iter >= a_projectileEnergy ) break;
    }

    if( index == 0 ) return( -1 );

    --index;

    double randomMax = a_random * m_levelsAndProbabilities.m_summedProbabilities[m_indices[index]];
    for( index = 0; index < m_levelsAndProbabilities.m_levels.size( ) - 1; ++index ) {
        if( m_levelsAndProbabilities.m_summedProbabilities[index] >= randomMax ) {
            break;
        }
    }
    return( m_levelsAndProbabilities.m_levels[index] );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void GRIN_inelasticForEnergy::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_VECTOR_INT( m_indices, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_thresholds, a_buffer, a_mode );
    m_levelsAndProbabilities.serialize( a_buffer, a_mode );
}

/*! \class GRIN_inelastic
 * This class represents GRIN inelastic continuum reaction data which has simulated levels.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_inelastic::GRIN_inelastic( ) {

}

/* *********************************************************************************************************//**
 * @param a_setupInfo                   [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param GRIN_continuumGammas          [in]    GIDI instance containing the GRIN capture data.
 ***********************************************************************************************************/

LUPI_HOST GRIN_inelastic::GRIN_inelastic( SetupInfo &a_setupInfo, GIDI::GRIN::GRIN_continuumGammas const &GRIN_continuumGammas ) :
        m_neutronIndex( a_setupInfo.m_neutronIndex ),
        m_neutronUserParticleIndex( -1 ),
        m_neutronMass( a_setupInfo.m_protare.projectileMass( ) ),

        m_targetIntid( a_setupInfo.m_protare.targetIntid( ) ),
        m_targetIndex( a_setupInfo.m_protare.targetIndex( ) ),
        m_targetUserParticleIndex( -1 ),
        m_targetMass( a_setupInfo.m_protare.targetMass( ) ) {

    PoPI::Database const &pops = GRIN_continuumGammas.pops( );
    GIDI::Suite const &inelasticIncidentEnergies = GRIN_continuumGammas.inelasticIncidentEnergies( );
    m_energies.reserve( inelasticIncidentEnergies.size( ) );
    m_inelasticForEnergy.reserve( inelasticIncidentEnergies.size( ) );
    for( std::size_t index = 0; index < inelasticIncidentEnergies.size( ); ++index ) {

        GIDI::GRIN::InelasticIncidentEnergy const *inelasticIncidentEnergy = inelasticIncidentEnergies.get<GIDI::GRIN::InelasticIncidentEnergy>( index );

        m_energies.push_back( inelasticIncidentEnergy->energy( ) );
        m_inelasticForEnergy.push_back( new GRIN_inelasticForEnergy( a_setupInfo, m_neutronMass, m_targetMass, pops, inelasticIncidentEnergy ) );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_inelastic::~GRIN_inelastic( ) {

    for( auto iter = m_inelasticForEnergy.begin( ); iter != m_inelasticForEnergy.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void GRIN_inelastic::setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {

    if( m_neutronIndex == a_particleIndex ) m_neutronUserParticleIndex = a_userParticleIndex;
    if( m_targetIndex == a_particleIndex ) m_targetUserParticleIndex = a_userParticleIndex;
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *
 * @param a_particleIntid       [in]    The PoPs intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void GRIN_inelastic::setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {

    if( PoPI::Intids::neutron == a_particleIntid ) m_neutronUserParticleIndex = a_userParticleIndex;
    if( m_targetIntid == a_particleIntid ) m_targetUserParticleIndex = a_userParticleIndex;
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void GRIN_inelastic::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_INT( m_neutronIndex, a_buffer, a_mode );
    DATA_MEMBER_INT( m_neutronUserParticleIndex, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_neutronMass, a_buffer, a_mode );

    DATA_MEMBER_INT( m_targetIntid, a_buffer, a_mode );
    DATA_MEMBER_INT( m_targetIndex, a_buffer, a_mode );
    DATA_MEMBER_INT( m_targetUserParticleIndex, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_targetMass, a_buffer, a_mode );

    DATA_MEMBER_VECTOR_DOUBLE( m_energies, a_buffer, a_mode );

    std::size_t vectorSize = m_inelasticForEnergy.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        m_inelasticForEnergy.resize( vectorSize, &a_buffer.m_placement );
        for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
            if( a_buffer.m_placement != nullptr ) {
                m_inelasticForEnergy[vectorIndex] = new(a_buffer.m_placement) GRIN_inelasticForEnergy;
                a_buffer.incrementPlacement( sizeof( GRIN_inelasticForEnergy ) ); }
            else {
               m_inelasticForEnergy[vectorIndex] = new GRIN_inelasticForEnergy;
            }
        } }
    else if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += m_inelasticForEnergy.internalSize( );
        a_buffer.incrementPlacement( sizeof( GRIN_inelasticForEnergy ) * vectorSize );
    }

    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_inelasticForEnergy[vectorIndex]->serialize( a_buffer, a_mode );
    }
}

/*! \class GRIN_captureToCompound
 * This class represents 
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_captureToCompound::GRIN_captureToCompound( ) {

}

/* *********************************************************************************************************//**
 * @param a_setupInfo           [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param a_compoundId          [in]    The GNDS PoPs' id of the compound level.
 ***********************************************************************************************************/

LUPI_HOST GRIN_captureToCompound::GRIN_captureToCompound( SetupInfo &a_setupInfo, PoPI::Database const &a_pops, std::string a_compoundId ) :
        m_index( a_setupInfo.m_stateNamesToIndices[a_compoundId] ),
        m_continuumIndices( ) {

        PoPI::Nuclide const &nuclide = a_pops.get<PoPI::Nuclide const>( a_compoundId );

        PoPI::GammaDecayData const &gammaDecayData = nuclide.gammaDecayData( );
        auto ids = gammaDecayData.ids( );
        auto probabilities1 = gammaDecayData.probabilities( );

        std::vector<int> levels;
        levels.reserve( ids.size( ) );
        std::vector<double> probabilities2;
        probabilities2.reserve( ids.size( ) );
        for( std::size_t index = 0; index < ids.size( ); ++index ) {
            PoPI::Nuclide const &daughter = a_pops.get<PoPI::Nuclide const>( ids[index] );
            if( daughter.kind( ) != PoPI_continuumChars ) continue;
            levels.push_back( a_setupInfo.m_stateNamesToIndices[ids[index]] );
            probabilities2.push_back( probabilities1[index] );
        }
        m_continuumIndices.set( levels, probabilities2 );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_captureToCompound::~GRIN_captureToCompound( ) {

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void GRIN_captureToCompound::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_INT( m_index, a_buffer, a_mode );
    m_continuumIndices.serialize( a_buffer, a_mode );
}

/*! \class GRIN_captureLevelProbability
 * This class represents 
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_captureLevelProbability::GRIN_captureLevelProbability( ) {

}

/* *********************************************************************************************************//**
 * @param a_setupInfo                   [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param a_captureLevelProbability     [in]    GIDI instance with the data.
 ***********************************************************************************************************/

LUPI_HOST GRIN_captureLevelProbability::GRIN_captureLevelProbability( SetupInfo &a_setupInfo, PoPI::Database const &a_pops,
                GIDI::GRIN::CaptureLevelProbability const *a_captureLevelProbability ) :
        m_knownLevelsAndProbabilities( a_setupInfo, a_pops, a_captureLevelProbability->table( ), false ) {

    auto capturePrimaryToContinua = LUPI::Misc::splitString( a_captureLevelProbability->capturePrimaryToContinua( ), true );
    m_captureToCompounds.reserve( capturePrimaryToContinua.size( ) );
    for( std::size_t i1 = 0; i1 < capturePrimaryToContinua.size( ); ++i1 ) {
        m_captureToCompounds.push_back( new GRIN_captureToCompound( a_setupInfo, a_pops, capturePrimaryToContinua[i1] ) );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_captureLevelProbability::~GRIN_captureLevelProbability( ) {

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void GRIN_captureLevelProbability::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    m_knownLevelsAndProbabilities.serialize( a_buffer, a_mode );

    std::size_t vectorSize = m_captureToCompounds.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        m_captureToCompounds.resize( vectorSize, &a_buffer.m_placement );
        for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
            if( a_buffer.m_placement != nullptr ) {
                m_captureToCompounds[vectorIndex] = new(a_buffer.m_placement) GRIN_captureToCompound;
                a_buffer.incrementPlacement( sizeof( GRIN_captureToCompound ) ); }
            else {
               m_captureToCompounds[vectorIndex] = new GRIN_captureToCompound;
            }
        } }
    else if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += m_captureToCompounds.internalSize( );
        a_buffer.incrementPlacement( sizeof( GRIN_captureToCompound ) * vectorSize );
    }

    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_captureToCompounds[vectorIndex]->serialize( a_buffer, a_mode );
    }
}

/*! \class GRIN_capture
 * Thiis class represents GRIN capture continuum reaction data which has simulated (i.e., modelled) levels.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_capture::GRIN_capture( ) :
        m_captureNeutronSeparationEnergy( 0.0 ),
        m_residualIntid( -1 ),
        m_residualIndex( -1 ),
        m_residualUserIndex( -1 ),
        m_residualMass( 0.0 ) {

}

/* *********************************************************************************************************//**
 * @param a_setupInfo                   [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param GRIN_continuumGammas          [in]    GIDI instance containing the GRIN capture data.
 ***********************************************************************************************************/

LUPI_HOST GRIN_capture::GRIN_capture( SetupInfo &a_setupInfo, GIDI::GRIN::GRIN_continuumGammas const &GRIN_continuumGammas ) :
        m_captureNeutronSeparationEnergy( GRIN_continuumGammas.captureNeutronSeparationEnergy( ).value( ) ),
        m_residualIntid( GRIN_continuumGammas.captureResidualIntid( ) ),
        m_residualIndex( GRIN_continuumGammas.captureResidualIndex( ) ),
        m_residualUserIndex( -1 ),
        m_residualMass( GRIN_continuumGammas.captureResidualMass( ) ) {

    PoPI::Database const &pops = GRIN_continuumGammas.pops( );
    GIDI::Suite const &captureLevelProbabilities = GRIN_continuumGammas.captureLevelProbabilities( );

    m_summedProbabilities.reserve( captureLevelProbabilities.size( ) );
    m_captureLevelProbabilities.reserve( captureLevelProbabilities.size( ) );
    double sum = 0.0;
    for( std::size_t index = 0; index < captureLevelProbabilities.size( ); ++index ) {
        GIDI::GRIN::CaptureLevelProbability const *captureLevelProbability = captureLevelProbabilities.get<GIDI::GRIN::CaptureLevelProbability const>( 0 );

        sum += captureLevelProbability->probabilty( );
        m_summedProbabilities.push_back( sum );
        m_captureLevelProbabilities.push_back( new GRIN_captureLevelProbability( a_setupInfo, pops, captureLevelProbability ) );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GRIN_capture::~GRIN_capture( ) {

    for( auto iter = m_captureLevelProbabilities.begin( ); iter != m_captureLevelProbabilities.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void GRIN_capture::setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {

    if( m_residualIndex == a_particleIndex ) m_residualUserIndex = a_userParticleIndex;
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *
 * @param a_particleIntid       [in]    The PoPs intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void GRIN_capture::setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {

    if( m_residualIntid == a_particleIntid ) m_residualUserIndex = a_userParticleIndex;
}


/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void GRIN_capture::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_DOUBLE( m_captureNeutronSeparationEnergy, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_summedProbabilities, a_buffer, a_mode );
    DATA_MEMBER_INT( m_residualIntid, a_buffer, a_mode );
    DATA_MEMBER_INT( m_residualIndex, a_buffer, a_mode );
    DATA_MEMBER_INT( m_residualUserIndex, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_residualMass, a_buffer, a_mode );

    std::size_t vectorSize = m_captureLevelProbabilities.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        m_captureLevelProbabilities.resize( vectorSize, &a_buffer.m_placement );
        for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
            if( a_buffer.m_placement != nullptr ) {
                m_captureLevelProbabilities[vectorIndex] = new(a_buffer.m_placement) GRIN_captureLevelProbability;
                a_buffer.incrementPlacement( sizeof( GRIN_captureLevelProbability ) ); }
            else {
               m_captureLevelProbabilities[vectorIndex] = new GRIN_captureLevelProbability;
            }
        } }
    else if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += m_captureLevelProbabilities.internalSize( );
        a_buffer.incrementPlacement( sizeof( GRIN_captureLevelProbability ) * vectorSize );
    }

    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_captureLevelProbabilities[vectorIndex]->serialize( a_buffer, a_mode );
    }
}

}
