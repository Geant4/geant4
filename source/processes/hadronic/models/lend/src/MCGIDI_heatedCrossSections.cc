/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "MCGIDI.hpp"

namespace MCGIDI {

static LUPI_HOST void checkZeroReaction( GIDI::Vector &vector, bool a_zeroReactions );
static LUPI_HOST GIDI::Vector collapseAndcheckZeroReaction( GIDI::Vector &a_vector, Transporting::MC const &a_settings, 
                GIDI::Transporting::Particles const &a_particles, double a_temperature, bool a_zeroReactions );
static void writeVector( FILE *a_file, std::string const &a_prefix, int a_offset, Vector<double> const &a_vector );

/*! \class HeatedReactionCrossSectionContinuousEnergy
 * Class to store a reaction's cross section.
 */

/* *********************************************************************************************************//**
 * Simple constructor needed for broadcasting.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE HeatedReactionCrossSectionContinuousEnergy::HeatedReactionCrossSectionContinuousEnergy( ) :
        m_offset( 0 ),
        m_threshold( 0.0 ),
        m_crossSections( ),
        m_URR_mode( Transporting::URR_mode::none ),
        m_URR_probabilityTables( nullptr ),
        m_ACE_URR_probabilityTables( nullptr ) {

}

/* *********************************************************************************************************//**
 * @param a_offset                  [in]    The offset of the first cross section point in the energy grid.
 * @param a_threshold               [in]    The threshold for the reaction.
 * @param a_crossSection            [in]    The cross section for the reaction.
 ***********************************************************************************************************/

LUPI_HOST HeatedReactionCrossSectionContinuousEnergy::HeatedReactionCrossSectionContinuousEnergy( int a_offset, double a_threshold, Vector<double> &a_crossSection ) :
        m_offset( a_offset ),
        m_threshold( a_threshold ),
        m_crossSections( a_crossSection.size( ) ),
        m_URR_mode( Transporting::URR_mode::none ),
        m_URR_probabilityTables( nullptr ),
        m_ACE_URR_probabilityTables( nullptr ) {

    int index = 0;              // This and next line needed as m_crossSections may be an instance of Vector<float>.
    for( auto iter = a_crossSection.begin( ); iter != a_crossSection.end( ); ++iter, ++index ) m_crossSections[index] = *iter;

}

/* *********************************************************************************************************//**
 * @param a_threshold                   [in]    The threshold for the reaction.
 * @param a_crossSection                [in]    The cross section for the reaction.
 * @param a_URR_probabilityTables       [in]    The pdf style URR probability tables.
 * @param a_ACE_URR_probabilityTables   [in]    The ACE style URR probability tables.
 ***********************************************************************************************************/

LUPI_HOST HeatedReactionCrossSectionContinuousEnergy::HeatedReactionCrossSectionContinuousEnergy( double a_threshold, 
                GIDI::Functions::Ys1d const &a_crossSection, Probabilities::ProbabilityBase2d *a_URR_probabilityTables,
                ACE_URR_probabilityTables *a_ACE_URR_probabilityTables ) :
        m_offset( a_crossSection.start( ) ),
        m_threshold( a_threshold ),
        m_crossSections( a_crossSection.Ys( ).size( ) ),
        m_URR_mode( Transporting::URR_mode::none ),
        m_URR_probabilityTables( a_URR_probabilityTables ),
        m_ACE_URR_probabilityTables( a_ACE_URR_probabilityTables ) {

    int index = 0;              // Next lines needed as m_crossSections may be an instance of Vector<float>.
    std::vector<double> const &Ys = a_crossSection.Ys( );
    for( auto iter = Ys.begin( ); iter != Ys.end( ); ++iter, ++index ) m_crossSections[index] = *iter;

    if( m_URR_probabilityTables != nullptr ) {
        m_URR_mode = Transporting::URR_mode::pdfs; }
    else if( m_ACE_URR_probabilityTables != nullptr ) {
        m_URR_mode = Transporting::URR_mode::ACE_URR_probabilityTables;
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE HeatedReactionCrossSectionContinuousEnergy::~HeatedReactionCrossSectionContinuousEnergy( ) {

    delete m_URR_probabilityTables;
    delete m_ACE_URR_probabilityTables;
}

/* *********************************************************************************************************//**
 * Returns the minimum energy for the unresolved resonance region (URR) domain. If no URR data present, returns -1.
 *
 * @return                              true is if *this* has a URR data.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedReactionCrossSectionContinuousEnergy::URR_domainMin( ) const {

    if( m_URR_probabilityTables != nullptr ) return( m_URR_probabilityTables->domainMin( ) );
    if( m_ACE_URR_probabilityTables != nullptr ) return( m_ACE_URR_probabilityTables->domainMin( ) );

    return( -1.0 );    
}

/* *********************************************************************************************************//**
 * Returns the maximum energy for the unresolved resonance region (URR) domain. If no URR data present, returns -1.
 *
 * @return                              true is if *this* has a URR data.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedReactionCrossSectionContinuousEnergy::URR_domainMax( ) const {

    if( m_URR_probabilityTables != nullptr ) return( m_URR_probabilityTables->domainMax( ) );
    if( m_ACE_URR_probabilityTables != nullptr ) return( m_ACE_URR_probabilityTables->domainMax( ) );

    return( -1.0 );    
}

/* *********************************************************************************************************//**
 * Returns the reactions cross section as a GIDI::Functions::XYs1d instance.
 *
 * @returns                     A GIDI::Functions::XYs1d instance.
 ***********************************************************************************************************/

LUPI_HOST GIDI::Functions::XYs1d HeatedReactionCrossSectionContinuousEnergy::crossSectionAsGIDI_XYs1d( double a_temperature,
                Vector<double> const &a_energies ) const {

    std::vector<double> energies = vectorToSTD_vector( a_energies );
    std::vector<double> crossSection( energies.size( ), 0.0 );
    for( std::size_t index = 0; index < m_crossSections.size( ); ++index ) crossSection[m_offset+index] = m_crossSections[index];

    return( GIDI::Functions::XYs1d( GIDI::Axes( ), ptwXY_interpolationLinLin, energies, crossSection, 0, a_temperature ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void HeatedReactionCrossSectionContinuousEnergy::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_INT( m_offset, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE(  m_threshold, a_buffer, a_mode  );
    DATA_MEMBER_VECTOR_FLOAT_OR_DOUBLE( m_crossSections, a_buffer, a_mode );
    m_URR_mode = serializeURR_mode( m_URR_mode,  a_buffer, a_mode );
    m_URR_probabilityTables = serializeProbability2d( a_buffer, a_mode, m_URR_probabilityTables );
    m_ACE_URR_probabilityTables = serializeACE_URR_probabilityTables( m_ACE_URR_probabilityTables, a_buffer, a_mode );
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

LUPI_HOST void HeatedReactionCrossSectionContinuousEnergy::print( LUPI_maybeUnused ProtareSingle const *a_protareSingle, std::string const &a_indent, 
                LUPI_maybeUnused std::string const &a_iFormat, LUPI_maybeUnused std::string const &a_energyFormat, std::string const &a_dFormat ) const {

    std::cout << a_indent << "# Offset = " << m_offset << std::endl;
    std::cout << a_indent << "# Threshold = " << m_threshold << std::endl;

    std::cout << a_indent << "# Number of cross section points = " << m_crossSections.size( ) << std::endl;
    for( auto iter = m_crossSections.begin( ); iter != m_crossSections.end( ); ++iter )
        std::cout << a_indent << LUPI::Misc::argumentsToString( a_dFormat.c_str( ), *iter ) << std::endl;
}

/*
============================================================
=================== ContinuousEnergyGain ===================
============================================================
*/

/* *********************************************************************************************************//**
 ***********************************************************************************************************/
LUPI_HOST_DEVICE ContinuousEnergyGain::ContinuousEnergyGain( ) :
        m_particleIntid( -1 ),
        m_particleIndex( -1 ),
        m_userParticleIndex( -1 ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST ContinuousEnergyGain::ContinuousEnergyGain( int a_particleIntid, int a_particleIndex, std::size_t a_size ) :
        m_particleIntid( a_particleIntid ),
        m_particleIndex( a_particleIndex ),
        m_userParticleIndex( -1 ),
        m_gain( a_size, static_cast<MCGIDI_FLOAT>( 0.0 ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST ContinuousEnergyGain &ContinuousEnergyGain::operator=( ContinuousEnergyGain const &a_continuousEnergyGain ) {

    m_particleIntid = a_continuousEnergyGain.particleIntid( );
    m_particleIndex = a_continuousEnergyGain.particleIndex( );
    m_userParticleIndex = a_continuousEnergyGain.userParticleIndex( );
    m_gain = a_continuousEnergyGain.gain( );

    return( *this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ContinuousEnergyGain::gain( int a_energy_index, double a_energy_fraction ) const {

    return( a_energy_fraction * m_gain[a_energy_index] + ( 1.0 - a_energy_fraction ) * m_gain[a_energy_index+1] );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ContinuousEnergyGain::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_INT( m_particleIntid, a_buffer, a_mode );
    DATA_MEMBER_INT( m_particleIndex, a_buffer, a_mode );
    DATA_MEMBER_INT( m_userParticleIndex, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_FLOAT_OR_DOUBLE( m_gain, a_buffer, a_mode );
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

LUPI_HOST void ContinuousEnergyGain::print( LUPI_maybeUnused ProtareSingle const *a_protareSingle, std::string const &a_indent, LUPI_maybeUnused std::string const &a_iFormat, 
                LUPI_maybeUnused std::string const &a_energyFormat, std::string const &a_dFormat ) const {

    std::cout << std::endl;
    std::cout << a_indent << "# Particle intid = " << m_particleIntid << std::endl;
    std::cout << a_indent << "# Particle index = " << m_particleIndex << std::endl;
    std::cout << a_indent << "# Use particle index = " << m_userParticleIndex << std::endl;

    std::cout << a_indent << "# Number of gains = " << m_gain.size( ) << std::endl;
    for( auto iter = m_gain.begin( ); iter != m_gain.end( ); ++iter )
        std::cout << a_indent << LUPI::Misc::argumentsToString( a_dFormat.c_str( ), *iter ) << std::endl;
}

/*
============================================================
=========== HeatedCrossSectionContinuousEnergy =============
============================================================
*/
LUPI_HOST_DEVICE HeatedCrossSectionContinuousEnergy::HeatedCrossSectionContinuousEnergy( ) :
        m_temperature( 0.0 ),
        m_hashIndices( ),
        m_energies( ),
        m_totalCrossSection( ),
        m_depositionEnergy( ),
        m_depositionMomentum( ),
        m_productionEnergy( ),
        m_gains( ),
        m_URR_mode( Transporting::URR_mode::none ),
        m_reactionsInURR_region( ),
        m_reactionCrossSections( ),
        m_ACE_URR_probabilityTables( nullptr ) {

}

/* *********************************************************************************************************//**
 * Fills in *this* with the requested temperature data.
 *
 * @param a_setupInfo                   [in]    Used internally when constructing a Protare to pass information to other components.
 * @param a_settings                    [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles                   [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 * @param a_domainHash                  [in]    The hash data used when looking up a cross section.
 * @param a_temperatureInfo             [in]    The list of temperatures to use.
 * @param a_reactions                   [in]    The list of reactions to use.
 * @param a_orphanProducts              [in]    The list of orphan products to use.
 * @param a_fixedGrid                   [in]    If true, the specified fixed grid is used; otherwise, grid in the file is used.
 * @param a_zeroReactions               [in]    Special case where no reaction in a protare is wanted so the first one is used but its cross section is set to 0.0 at all energies.
 ***********************************************************************************************************/

LUPI_HOST HeatedCrossSectionContinuousEnergy::HeatedCrossSectionContinuousEnergy( SetupInfo &a_setupInfo, Transporting::MC const &a_settings, 
                GIDI::Transporting::Particles const &a_particles, DomainHash const &a_domainHash, GIDI::Styles::TemperatureInfo const &a_temperatureInfo, 
                std::vector<GIDI::Reaction const *> const &a_reactions, std::vector<GIDI::Reaction const *> const &a_orphanProducts, bool a_fixedGrid,
                bool a_zeroReactions ) :
        m_temperature( a_temperatureInfo.temperature( ).value( ) ),
        m_hashIndices( ),
        m_energies( ),
        m_totalCrossSection( ),
        m_depositionEnergy( ),
        m_depositionMomentum( ),
        m_productionEnergy( ),
        m_gains( ),
        m_URR_mode( Transporting::URR_mode::none ),
        m_reactionsInURR_region( ),
        m_reactionCrossSections( ),
        m_ACE_URR_probabilityTables( nullptr ) {

    bool isPhotoAtomic = a_setupInfo.m_protare.isPhotoAtomic( );
    std::string label( a_temperatureInfo.griddedCrossSection( ) );
    std::string URR_label( a_temperatureInfo.URR_probabilityTables( ) );

    GIDI::Styles::GriddedCrossSection const &griddedCrossSectionStyle = 
            static_cast<GIDI::Styles::GriddedCrossSection const &>( *a_settings.styles( )->get<GIDI::Styles::Base>( label ) );
    GIDI::Grid const &grid = griddedCrossSectionStyle.grid( );

    std::vector<double> const *energiesPointer;
    std::vector<double> const &energies = grid.data( ).vector();
    std::vector<double> const &fixedGridPoints = a_settings.fixedGridPoints( );
    std::vector<int> fixedGridIndices( fixedGridPoints.size( ) );
    if( a_fixedGrid ) {
        for( std::size_t i1 = 0; i1 < fixedGridPoints.size( ); ++i1 ) {
            fixedGridIndices[i1] = binarySearchVector( fixedGridPoints[i1], energies );
        }
        energiesPointer = &fixedGridPoints;
        m_energies = fixedGridPoints; }
    else {
        energiesPointer = &energies;
        m_energies = energies;
    }

    m_hashIndices = a_domainHash.map( m_energies );

    int reactionIndex = 0;
    GIDI::Axes axes;
    std::vector<double> dummy;
    GIDI::Functions::Ys1d totalCrossSection( axes, ptwXY_interpolationLinLin, 0, dummy );
    GIDI::Functions::Ys1d fixedGridCrossSection( axes, ptwXY_interpolationLinLin, 0, fixedGridPoints );
    m_reactionCrossSections.resize( a_reactions.size( ) );
    m_URR_mode = Transporting::URR_mode::none;
    for( std::vector<GIDI::Reaction const *>::const_iterator reactionIter = a_reactions.begin( ); reactionIter != a_reactions.end( ); ++reactionIter, ++reactionIndex ) {
        GIDI::Suite const &reactionCrossSectionSuite = (*reactionIter)->crossSection( );
        GIDI::Functions::Ys1d const *reactionCrossSection3 = reactionCrossSectionSuite.get<GIDI::Functions::Ys1d>( label );
        GIDI::Functions::Ys1d *reactionCrossSectionZeroReactions = nullptr;
        if( a_zeroReactions ) {
            reactionCrossSectionZeroReactions = new GIDI::Functions::Ys1d( *reactionCrossSection3 );
            for( std::size_t index = 0; index < reactionCrossSectionZeroReactions->size( ); ++index ) reactionCrossSectionZeroReactions->set( index, 0.0 );
            reactionCrossSection3 = reactionCrossSectionZeroReactions;
        }

        Probabilities::ProbabilityBase2d *URR_probabilityTables = nullptr;
        if( ( a_settings._URR_mode( ) == Transporting::URR_mode::pdfs ) && ( URR_label != "" ) ) {
            if( reactionCrossSectionSuite.has( URR_label ) ) {
                GIDI::Functions::URR_probabilityTables1d const &URR_probability_tables1d( *reactionCrossSectionSuite.get<GIDI::Functions::URR_probabilityTables1d>( URR_label ) );
                URR_probabilityTables = Probabilities::parseProbability2d( URR_probability_tables1d.function2d( ), nullptr );
                m_URR_mode = Transporting::URR_mode::pdfs;
            }
        }

        ACE_URR_probabilityTables *ACE_URR_probabilityTables1 = nullptr;
        if( a_settings._URR_mode( ) == Transporting::URR_mode::ACE_URR_probabilityTables ) {
            auto URR_iter = a_setupInfo.m_ACE_URR_probabilityTablesFromGIDI.find( URR_label );
            if( URR_iter != a_setupInfo.m_ACE_URR_probabilityTablesFromGIDI.end( ) ) {
                auto ACE_URR_probabilityTablesIter = (*URR_iter).second->m_ACE_URR_probabilityTables.find( (*reactionIter)->label( ) );
                if( ACE_URR_probabilityTablesIter != (*URR_iter).second->m_ACE_URR_probabilityTables.end( ) ) {
                    ACE_URR_probabilityTables1 = (*ACE_URR_probabilityTablesIter).second;
                    (*ACE_URR_probabilityTablesIter).second = nullptr;          // Set to nullptr so destructor of m_ACE_URR_probabilityTables does not delete.
                    m_URR_mode = Transporting::URR_mode::ACE_URR_probabilityTables;
                }
            }
        }

        if( a_fixedGrid ) {
            GIDI::Functions::Ys1d *reactionCrossSection4 = &fixedGridCrossSection;

            int start = 0;

            if( energies[reactionCrossSection3->start( )] > fixedGridPoints[0] ) {
                start = binarySearchVector( energies[reactionCrossSection3->start( )], fixedGridPoints ) + 1;
            }

            for( int i1 = 0; i1 < start; ++i1 ) reactionCrossSection4->set( i1, 0.0 );
            for( int i1 = start; i1 < static_cast<int>( fixedGridPoints.size( ) ); ++i1 ) {
                int index = fixedGridIndices[i1];
                double fraction = ( fixedGridPoints[i1] - energies[index] ) / ( energies[index+1] - energies[index] );

                index -= reactionCrossSection3->start( );
                reactionCrossSection4->set( i1, ( 1.0 - fraction ) * (*reactionCrossSection3)[index] + fraction * (*reactionCrossSection3)[index+1] );
            }
            reactionCrossSection3 = &fixedGridCrossSection;
        }

        m_reactionCrossSections[reactionIndex] = new HeatedReactionCrossSectionContinuousEnergy( (*reactionIter)->crossSectionThreshold( ), 
                *reactionCrossSection3, URR_probabilityTables, ACE_URR_probabilityTables1 );
        totalCrossSection += *reactionCrossSection3;

        delete reactionCrossSectionZeroReactions;
    }
    m_totalCrossSection.resize( totalCrossSection.length( ), 0.0 );
    for( std::size_t i1 = 0; i1 < totalCrossSection.size( ); ++i1 ) m_totalCrossSection[i1+totalCrossSection.start()] = totalCrossSection[i1];

    if( hasURR_probabilityTables( ) ) {
        std::vector<int> reactions_in_URR_region;

        for( reactionIndex = 0; reactionIndex < numberOfReactions( ); ++reactionIndex ) {
            if( m_reactionCrossSections[reactionIndex]->threshold( ) < URR_domainMax( ) ) {
                reactions_in_URR_region.push_back( reactionIndex );
            }
        }

        m_reactionsInURR_region.resize( reactions_in_URR_region.size( ) );
        for( std::size_t i1 = 0; i1 < reactions_in_URR_region.size( ); ++i1 ) m_reactionsInURR_region[i1] = reactions_in_URR_region[i1];

        if( a_settings._URR_mode( ) == Transporting::URR_mode::ACE_URR_probabilityTables ) {
            auto URR_iter = a_setupInfo.m_ACE_URR_probabilityTablesFromGIDI.find( URR_label );
            if( URR_iter != a_setupInfo.m_ACE_URR_probabilityTablesFromGIDI.end( ) ) {
                auto ACE_URR_probabilityTablesIter = (*URR_iter).second->m_ACE_URR_probabilityTables.find( "total" );
                if( ACE_URR_probabilityTablesIter != (*URR_iter).second->m_ACE_URR_probabilityTables.end( ) ) {
                    m_ACE_URR_probabilityTables = (*ACE_URR_probabilityTablesIter).second;
                    (*ACE_URR_probabilityTablesIter).second = nullptr;          // Set to nullptr so destructor of m_ACE_URR_probabilityTables does not delete.
                }
            }
        }
    }

    if( a_settings.addExpectedValueData( ) ) {
        m_depositionEnergy.resize( totalCrossSection.length( ), 0.0 );
        m_depositionMomentum.resize( totalCrossSection.length( ), 0.0 );
        m_productionEnergy.resize( totalCrossSection.length( ), 0.0 );

        m_gains.resize( a_particles.particles( ).size( ) );
        int i2 = 0;
        int projectileGainIndex = -1;
        int photonGainIndex = -1;
        for( std::map<std::string, GIDI::Transporting::Particle>::const_iterator particle = a_particles.particles( ).begin( ); particle != a_particles.particles( ).end( );
                        ++particle, ++i2 ) {
            int particleIntid = a_setupInfo.m_particleIntids[particle->first];
            int particleIndex = a_setupInfo.m_particleIndices[particle->first];

            if( particleIntid == a_setupInfo.m_protare.projectileIntid( ) ) projectileGainIndex = i2;
            m_gains[i2] = new ContinuousEnergyGain( particleIntid, particleIndex, totalCrossSection.length( ) );
            if( particle->first == PoPI::IDs::photon ) photonGainIndex = i2;
        }

        std::vector< std::vector<double> > gains( a_particles.particles( ).size( ) );
        for( std::size_t reactionIndex2 = 0; reactionIndex2 < a_reactions.size( ) + a_orphanProducts.size( ); ++reactionIndex2 ) {

            int offset = 0;
            std::vector<double> deposition_energy( m_energies.size( ), 0.0 );
            std::vector<double> deposition_momentum( m_energies.size( ), 0.0 );
            std::vector<double> production_energy( m_energies.size( ), 0.0 );
            for( std::size_t i1 = 0; i1 < gains.size( ); ++i1 ) gains[i1] = std::vector<double>( m_energies.size( ), 0.0 );

            HeatedReactionCrossSectionContinuousEnergy *MCGIDI_reaction_cross_section = nullptr;
            GIDI::Reaction const *reaction = nullptr;
            GIDI::Functions::Ys1d const *reactionCrossSection = nullptr;                // If nullptr, reaction is an orphanProduct node.
            GIDI::Functions::Function1dForm const *available_energy = nullptr;
            GIDI::Functions::Function1dForm const *available_momentum = nullptr;
            if( reactionIndex2 < a_reactions.size( ) ) {
                reaction = a_reactions[reactionIndex2];
                available_energy = reaction->availableEnergy( ).get<GIDI::Functions::Function1dForm>( 0 );
                available_momentum = reaction->availableMomentum( ).get<GIDI::Functions::Function1dForm>( 0 );
                MCGIDI_reaction_cross_section = m_reactionCrossSections[reactionIndex2];
                offset = MCGIDI_reaction_cross_section->offset( ); }
            else {
                reaction = a_orphanProducts[reactionIndex2-a_reactions.size( )];
                GIDI::Suite const &reactionCrossSectionSuite = reaction->crossSection( );
                reactionCrossSection = reactionCrossSectionSuite.get<GIDI::Functions::Ys1d>( label );
                offset = reactionCrossSection->start( );
            }
            if( a_settings.useSlowerContinuousEnergyConversion( ) ) {                   // Old way which is slow as it does one energy at a time.
                for( std::size_t energy_index = (std::size_t) offset; energy_index < m_energies.size( ); ++energy_index ) {
                    double energy = m_energies[energy_index];

                    if( reactionCrossSection == nullptr ) {
                        if( isPhotoAtomic ) {                                           // Treat as Q = 0.0 since 2 photons will be emitted.
                            deposition_energy[energy_index] = energy; }
                        else {
                            deposition_energy[energy_index] = available_energy->evaluate( energy );

                            double Q = deposition_energy[energy_index] - energy;        // Should use Q node to get this.
                            if( fabs( Q ) < 1e-12 * deposition_energy[energy_index] )
                                    Q = 0.0;                                            // Probably 0.0 due to rounding errors.
                            production_energy[energy_index] = Q;
                        }
                        deposition_momentum[energy_index] = available_momentum->evaluate( energy );
                    }

                    int i1 = 0;
                    for( std::map<std::string, GIDI::Transporting::Particle>::const_iterator particle = a_particles.particles( ).begin( ); particle != a_particles.particles( ).end( ); 
                            ++particle, ++i1 ) {
                        double product_energy, product_momentum, product_gain;

                        if( particle->first == PoPI::IDs::electron ) continue;      // As of this coding, electrons are not complete in GNDS files.
                                                                                    // When they are, this statement can be removed.
                        if( ( reactionCrossSection != nullptr ) && ( particle->first != PoPI::IDs::photon ) ) continue;

                        if( reaction->isPairProduction( ) && ( particle->first == PoPI::IDs::photon ) ) {
                            product_energy = 2.0 * PoPI_electronMass_MeV_c2;        // Assumes energy unit is MeV.
                            product_momentum = 0.0;
                            product_gain = 2.0; }
                        else {
                            reaction->continuousEnergyProductData( a_settings, particle->first, energy, product_energy, product_momentum, 
                                    product_gain, true );
                        }
                        if( i1 == projectileGainIndex ) --product_gain;

                        deposition_energy[energy_index] -= product_energy;
                        deposition_momentum[energy_index] -= product_momentum;
                        gains[i1][energy_index] = product_gain;
                    }
                } }
            else {                                                                  // New way which is hopefully faster.
                if( reactionCrossSection == nullptr ) {
                    if( isPhotoAtomic ) {                                           // Treat as Q = 0.0 since 2 photons will be emitted.
                        for( std::size_t energyIndex = (std::size_t) offset; energyIndex < m_energies.size( ); ++energyIndex ) {
                            deposition_energy[energyIndex] = m_energies[energyIndex];
                        } }
                    else {
                        available_energy->mapToXsAndAdd( offset, *energiesPointer, deposition_energy, 1.0 );
                        for( std::size_t energyIndex = (std::size_t) offset; energyIndex < m_energies.size( ); ++energyIndex ) {
                            double Q = deposition_energy[energyIndex] - m_energies[energyIndex];
                            if( fabs( Q ) < 1e-12 * deposition_energy[energyIndex] )
                                    Q = 0.0;                                            // Probably 0.0 due to rounding errors.
                            production_energy[energyIndex] = Q;
                        }
                    }
                    available_momentum->mapToXsAndAdd( offset, *energiesPointer, deposition_momentum, 1.0 );
                }

                int i1 = 0;
                for( std::map<std::string, GIDI::Transporting::Particle>::const_iterator particle = a_particles.particles( ).begin( );
                            particle != a_particles.particles( ).end( ); ++particle, ++i1 ) {
                    if( particle->first == PoPI::IDs::electron ) continue;      // As of this coding, electrons are not complete in GNDS files.
                                                                                // When they are, this statement can be removed.
                    if( ( reactionCrossSection != nullptr ) && ( particle->first != PoPI::IDs::photon ) ) continue;

                    if( reaction->isPairProduction( ) && ( particle->first == PoPI::IDs::photon ) ) {
                        for( std::size_t energyIndex = (std::size_t) offset; energyIndex < m_energies.size( ); ++energyIndex ) {
                            deposition_energy[energyIndex] -= 2.0 * PoPI_electronMass_MeV_c2;    // Assumes energy unit is MeV.
                            gains[i1][energyIndex] = 2.0;
                        } }
                    else {
                        reaction->mapContinuousEnergyProductData( a_settings, particle->first, *energiesPointer, offset, deposition_energy, 
                                deposition_momentum, gains[i1], true );
                    }

                    if( i1 == projectileGainIndex ) {
                        for( std::size_t energyIndex = (std::size_t) offset; energyIndex < m_energies.size( ); ++energyIndex ) --gains[i1][energyIndex];
                    }
                }
            }

            if( a_particles.hasParticle( PoPI::IDs::photon ) ) {
                if( a_setupInfo.m_initialStateIndices.find( reaction->label( ) ) != a_setupInfo.m_initialStateIndices.end( ) ) {
                    int initialStateIndex = a_setupInfo.m_initialStateIndices[reaction->label( )];
                    if( initialStateIndex >= 0 ) {
                        NuclideGammaBranchStateInfo *nuclideGammaBranchStateInfo = a_setupInfo.m_protare.nuclideGammaBranchStateInfos( )[initialStateIndex];
                        double multiplicity = nuclideGammaBranchStateInfo->multiplicity( );
                        double averageGammaEnergy = nuclideGammaBranchStateInfo->averageGammaEnergy( );

                        for( std::size_t energyIndex = (std::size_t) offset; energyIndex < m_energies.size( ); ++energyIndex ) {
                            deposition_energy[energyIndex] -= averageGammaEnergy;
                            if( photonGainIndex >= 0 ) gains[photonGainIndex][energyIndex] += multiplicity;
                        }
                    }
                }
            }

            double crossSection = 0.0;
            for( std::size_t energyIndex = (std::size_t) offset; energyIndex < m_energies.size( ); ++energyIndex ) {
                if( reactionCrossSection == nullptr ) {
                    crossSection = MCGIDI_reaction_cross_section->crossSection( energyIndex ); }
                else {
                    crossSection = (*reactionCrossSection)[energyIndex-offset];
                }

                m_depositionEnergy[energyIndex] += crossSection * deposition_energy[energyIndex];
                m_depositionMomentum[energyIndex] += crossSection * deposition_momentum[energyIndex];
                m_productionEnergy[energyIndex] += crossSection * production_energy[energyIndex];
                for( std::size_t i1 = 0; i1 < m_gains.size( ); ++i1 ) {
                    m_gains[i1]->adjustGain( energyIndex, crossSection * gains[i1][energyIndex] );
                }
            }
        }
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE HeatedCrossSectionContinuousEnergy::~HeatedCrossSectionContinuousEnergy( ) {

    for( auto iter = m_reactionCrossSections.begin( ); iter < m_reactionCrossSections.end( ); ++iter ) delete *iter;
    for( auto iter = m_gains.begin( ); iter < m_gains.end( ); ++iter ) delete *iter;
    delete m_ACE_URR_probabilityTables;
}


/* *********************************************************************************************************//**
 * This function returns the index in *a_energies* where *a_energy* lies between the returned index and the next index.
 * The returned index must lie between a_hashIndices[a_hashIndex] and a_hashIndices[a_hashIndex+1].
 * If *a_energy* is below the domain of *a_energies*, 0 is returned. If *a_energy* is above the domain of *a_energies*,
 * the size of *a_energies* minus 2 is returned.
 * The argument *a_energyFraction* the weight for the energy at the returned index with the next index getting weighting 1 minus
 * *a_energyFraction*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_energy              [in]    The energy whose index is requested.
 * @param a_energyFraction      [in]    This represents the weighting to apply to the two bounding energies.
 *
 * @return                              The index bounding *a_energy* in the member *m_energies*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE int HeatedCrossSectionContinuousEnergy::evaluationInfo( int a_hashIndex, double a_energy, double *a_energyFraction ) const {

    return Sampling::evaluationForHashIndex( a_hashIndex, m_hashIndices, a_energy, m_energies, a_energyFraction );
}

/* *********************************************************************************************************//**
 * Returns true if *this* has a unresolved resonance region (URR) data and false otherwise.
 *
 * @return                              true is if *this* has a URR data.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE bool HeatedCrossSectionContinuousEnergy::hasURR_probabilityTables( ) const {

    return( m_URR_mode != Transporting::URR_mode::none );
}

/* *********************************************************************************************************//**
 * Returns the minimum energy for the unresolved resonance region (URR) domain. If no URR data present, returns -1.
 *
 * @return                              true is if *this* has a URR data.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::URR_domainMin( ) const {

    if( m_ACE_URR_probabilityTables != nullptr ) return( m_ACE_URR_probabilityTables->domainMin( ) );

    for( std::size_t i1 = 0; i1 < m_reactionCrossSections.size( ); ++i1 ) {
        HeatedReactionCrossSectionContinuousEnergy *reactionCrossSection = m_reactionCrossSections[i1];

        if( reactionCrossSection->hasURR_probabilityTables( ) ) return( reactionCrossSection->URR_domainMin( ) );
    }

    return( -1.0 );
}

/* *********************************************************************************************************//**
 * Returns the maximum energy for the unresolved resonance region (URR) domain. If no URR data present, returns -1.
 *
 * @return                              true is if *this* has a URR data.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::URR_domainMax( ) const {

    if( m_ACE_URR_probabilityTables != nullptr ) return( m_ACE_URR_probabilityTables->domainMax( ) );

    for( std::size_t i1 = 0; i1 < m_reactionCrossSections.size( ); ++i1 ) {
        HeatedReactionCrossSectionContinuousEnergy *reactionCrossSection = m_reactionCrossSections[i1];

        if( reactionCrossSection->hasURR_probabilityTables( ) ) return( reactionCrossSection->URR_domainMax( ) );
    }

    return( -1.0 );
}
/*
=========================================================
*/
LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::crossSection( URR_protareInfos const &a_URR_protareInfos, int a_URR_index, 
                int a_hashIndex, double a_energy, LUPI_maybeUnused bool a_sampling ) const {

    double energy_fraction;
    int energy_index = evaluationInfo( a_hashIndex, a_energy, &energy_fraction );

    if( a_URR_index >= 0 ) {
        URR_protareInfo const &URR_protare_info = a_URR_protareInfos[a_URR_index];

        if( URR_protare_info.m_inURR ) {
            double cross_section = 0.0;
            for( std::size_t i1 = 0; i1 < m_reactionsInURR_region.size( ); ++i1 ) {
                cross_section += reactionCrossSection2( m_reactionsInURR_region[i1], a_URR_protareInfos, a_URR_index, a_energy, 
                        energy_index, energy_fraction, false );
            }
            return( cross_section );
        }
    }

    return( energy_fraction * m_totalCrossSection[energy_index] + ( 1.0 - energy_fraction ) * m_totalCrossSection[energy_index+1] );
}
/*
=========================================================
*/
LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, 
                int a_URR_index, int a_hashIndex, double a_energy, LUPI_maybeUnused bool a_sampling ) const {

    double energyFraction;

    int energyIndex = evaluationInfo( a_hashIndex, a_energy, &energyFraction );
    return( reactionCrossSection2( a_reactionIndex, a_URR_protareInfos, a_URR_index, a_energy, energyIndex, energyFraction ) );
}
/*
=========================================================
*/
LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::reactionCrossSection2( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, 
                int a_URR_index, double a_energy, int a_energyIndex, double a_energyFraction, LUPI_maybeUnused bool a_sampling ) const {

    HeatedReactionCrossSectionContinuousEnergy const &reaction = *m_reactionCrossSections[a_reactionIndex];
    double URR_cross_section_factor = 1.0;

    if( a_URR_index >= 0 ) {
        URR_protareInfo const &URR_protare_info = a_URR_protareInfos[a_URR_index];
        if( URR_protare_info.m_inURR ) {
            if( m_URR_mode == Transporting::URR_mode::pdfs ) {
                if( reaction.URR_probabilityTables( ) != nullptr )
                    URR_cross_section_factor = reaction.URR_probabilityTables( )->sample( a_energy, URR_protare_info.m_rng_Value, []() -> double { return 0.0; } ); }
            else if( m_URR_mode == Transporting::URR_mode::ACE_URR_probabilityTables ) {
                if( reaction._ACE_URR_probabilityTables( ) != nullptr ) 
                    URR_cross_section_factor = reaction._ACE_URR_probabilityTables( )->sample( a_energy, URR_protare_info.m_rng_Value );
            }
        }
    }

    return( URR_cross_section_factor * ( a_energyFraction * reaction.crossSection( a_energyIndex ) + ( 1.0 - a_energyFraction ) * reaction.crossSection( a_energyIndex + 1 ) ) );
}

/*
=========================================================
*/
LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_URR_index, 
                double a_energy_in ) const {

    int energyIndex = binarySearchVector( a_energy_in, m_energies );
    double energyFraction;

    if( energyIndex < 0 ) {
        if( energyIndex == -1 ) {
            energyIndex = static_cast<int>( m_energies.size( ) ) - 2;
            energyFraction = 0.0; }
        else {
            energyIndex = 0;
            energyFraction = 1.0;
        } }
    else {
        energyFraction = ( m_energies[energyIndex+1] - a_energy_in ) / ( m_energies[energyIndex+1] - m_energies[energyIndex] );
    }

    return( reactionCrossSection2( a_reactionIndex, a_URR_protareInfos, a_URR_index, a_energy_in, energyIndex, energyFraction, false ) );
}

/* *********************************************************************************************************//**
 * Returns the total cross section as a GIDI::Functions::XYs1d instance.
 *
 * @returns                     A GIDI::Functions::XYs1d instance.
 ***********************************************************************************************************/

LUPI_HOST GIDI::Functions::XYs1d HeatedCrossSectionContinuousEnergy::crossSectionAsGIDI_XYs1d( ) const {

    std::vector<double> energies = vectorToSTD_vector( m_energies );
    std::vector<double> crossSection = vectorToSTD_vector( m_totalCrossSection );

    return( GIDI::Functions::XYs1d( GIDI::Axes( ), ptwXY_interpolationLinLin, energies, crossSection, 0, m_temperature ) );
}

/* *********************************************************************************************************//**
 * Returns the reaction's cross section as GIDI::Functions::XYs1d instance.
 *
 * @param a_reactionIndex       [in]    Specifies the indexs of the reaction whose cross section is requested.
 *
 * @returns                             A GIDI::Functions::XYs1d instance.
 ***********************************************************************************************************/

LUPI_HOST GIDI::Functions::XYs1d HeatedCrossSectionContinuousEnergy::reactionCrossSectionAsGIDI_XYs1d( int a_reactionIndex ) const {

    return( m_reactionCrossSections[a_reactionIndex]->crossSectionAsGIDI_XYs1d( m_temperature, m_energies ) );
}

/* *********************************************************************************************************//**
 * Returns the index of a sampled reaction for a projectile with energy *a_energy* and total cross section
 * *a_crossSection*. Random numbers are obtained via *a_userrng* and *a_rngState*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_energy              [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::depositionEnergy( int a_hashIndex, double a_energy ) const {

    double energy_fraction;
    int energy_index = evaluationInfo( a_hashIndex, a_energy, &energy_fraction );

    return( energy_fraction * m_depositionEnergy[energy_index] + ( 1.0 - energy_fraction ) * m_depositionEnergy[energy_index+1] );
}

/* *********************************************************************************************************//**
 * Returns the index of a sampled reaction for a projectile with energy *a_energy* and total cross section
 * *a_crossSection*. Random numbers are obtained via *a_userrng* and *a_rngState*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_energy              [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::depositionMomentum( int a_hashIndex, double a_energy ) const {

    double energy_fraction;
    int energy_index = evaluationInfo( a_hashIndex, a_energy, &energy_fraction );

    return( energy_fraction * m_depositionMomentum[energy_index] + ( 1.0 - energy_fraction ) * m_depositionMomentum[energy_index+1] );
}

/* *********************************************************************************************************//**
 * Returns the index of a sampled reaction for a projectile with energy *a_energy* and total cross section
 * *a_crossSection*. Random numbers are obtained via *a_userrng* and *a_rngState*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_energy              [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::productionEnergy( int a_hashIndex, double a_energy ) const {

    double energy_fraction;
    int energy_index = evaluationInfo( a_hashIndex, a_energy, &energy_fraction );

    return( energy_fraction * m_productionEnergy[energy_index] + ( 1.0 - energy_fraction ) * m_productionEnergy[energy_index+1] );
}

/* *********************************************************************************************************//**
 * Returns the gain for the particle with index *a_particleIndex* for projectile energy *a_energy*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_particleIndex       [in]    The index of the particle whose gain is requested.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::gain( int a_hashIndex, double a_energy, int a_particleIndex ) const {

    double energy_fraction;
    int energy_index = evaluationInfo( a_hashIndex, a_energy, &energy_fraction );

    for( std::size_t i1 = 0; i1 < m_gains.size( ); ++i1 ) {
        if( a_particleIndex == m_gains[i1]->particleIndex( ) ) return( m_gains[i1]->gain( energy_index, energy_fraction ) );
    }

    return( 0.0 );
}

/* *********************************************************************************************************//**
 * Returns the gain for the particle with intid *a_particleIntid* for projectile energy *a_energy*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_particleIntid       [in]    The intid of the particle whose gain is requested.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionContinuousEnergy::gainViaIntid( int a_hashIndex, double a_energy, int a_particleIntid ) const {

    double energy_fraction;
    int energy_index = evaluationInfo( a_hashIndex, a_energy, &energy_fraction );

    for( std::size_t i1 = 0; i1 < m_gains.size( ); ++i1 ) {
        if( a_particleIntid == m_gains[i1]->particleIntid( ) ) return( m_gains[i1]->gain( energy_index, energy_fraction ) );
    }

    return( 0.0 );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionContinuousEnergy::setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {

    for( auto iter = m_gains.begin( ); iter != m_gains.end( ); ++iter ) (*iter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *
 * @param a_particleIntid       [in]    The intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionContinuousEnergy::setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {

    for( auto iter = m_gains.begin( ); iter != m_gains.end( ); ++iter ) (*iter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void HeatedCrossSectionContinuousEnergy::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_DOUBLE( m_temperature, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_hashIndices, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_energies, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_FLOAT_OR_DOUBLE( m_totalCrossSection, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_FLOAT_OR_DOUBLE( m_depositionEnergy, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_FLOAT_OR_DOUBLE( m_depositionMomentum, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_FLOAT_OR_DOUBLE( m_productionEnergy, a_buffer, a_mode );
    m_URR_mode = serializeURR_mode( m_URR_mode,  a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_reactionsInURR_region, a_buffer, a_mode );
    m_ACE_URR_probabilityTables = serializeACE_URR_probabilityTables( m_ACE_URR_probabilityTables, a_buffer, a_mode );

    std::size_t vectorSize = m_reactionCrossSections.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_reactionCrossSections.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_reactionCrossSections.internalSize();
    for( std::size_t memberIndex = 0; memberIndex < vectorSize; ++memberIndex ) {
        if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
            if( a_buffer.m_placement != nullptr ) {
                m_reactionCrossSections[memberIndex] = new(a_buffer.m_placement) HeatedReactionCrossSectionContinuousEnergy;
                a_buffer.incrementPlacement( sizeof( HeatedReactionCrossSectionContinuousEnergy ) ); }
            else {
                m_reactionCrossSections[memberIndex] = new HeatedReactionCrossSectionContinuousEnergy;
            }
        }
        if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof( HeatedReactionCrossSectionContinuousEnergy ) );
        }
        m_reactionCrossSections[memberIndex]->serialize( a_buffer, a_mode );
    }

    vectorSize = m_gains.size( );
    vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_gains.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_gains.internalSize();
    for( std::size_t memberIndex = 0; memberIndex < vectorSize; ++memberIndex ) {
        if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
            if( a_buffer.m_placement != nullptr ) {
                m_gains[memberIndex] = new(a_buffer.m_placement) ContinuousEnergyGain;
                a_buffer.incrementPlacement( sizeof( ContinuousEnergyGain ) ); }
            else {
                m_gains[memberIndex] = new ContinuousEnergyGain;
            }
        }

        if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof( ContinuousEnergyGain ) );
        }
        m_gains[memberIndex]->serialize( a_buffer, a_mode );
    }
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

LUPI_HOST void HeatedCrossSectionContinuousEnergy::print( ProtareSingle const *a_protareSingle, std::string const &a_indent, 
                std::string const &a_iFormat, std::string const &a_energyFormat, std::string const &a_dFormat ) const {

    char const *dFormat = a_dFormat.c_str( );
    std::string indent2 = a_indent + "  ";
    std::string lineFormat = indent2 + a_energyFormat + " " + a_dFormat + " " + a_dFormat + " " + a_dFormat + " " + a_dFormat;

    std::cout << std::endl;
    std::cout << a_indent << "# Temperature = " << LUPI::Misc::doubleToString3( dFormat, m_temperature, true ) << std::endl;

    std::cout << a_indent << "# Number of hash indices = " << m_hashIndices.size( ) << std::endl;
    for( auto iter = m_hashIndices.begin( ); iter != m_hashIndices.end( ); ++iter )
        std::cout << indent2 << LUPI::Misc::argumentsToString( a_iFormat.c_str( ), *iter ) << std::endl;
    std::cout << std::endl;
    std::size_t energySize = m_energies.size( );
    std::cout << a_indent << "# Number of energies = " << m_energies.size( ) << std::endl;
    std::cout << "#      projectile              total              deposition           deposition          production" << std::endl;
    std::cout << "#        energy             cross section           energy              momentum             energy"<< std::endl;
    std::cout << "#====================================================================================================="<< std::endl;
    for( std::size_t index = 0; index != energySize; ++index ) {
        std::cout << LUPI::Misc::argumentsToString( lineFormat.c_str( ), m_energies[index], m_totalCrossSection[index], m_depositionEnergy[index], 
                m_depositionMomentum[index], m_productionEnergy[index] ) << std::endl;
    }

    for( auto iter = m_gains.begin( ); iter != m_gains.end( ); ++iter ) (*iter)->print( a_protareSingle, a_indent, a_iFormat, a_energyFormat, a_dFormat );

    std::cout << std::endl;
    std::cout << "# URR mode is ";
    if( m_URR_mode == Transporting::URR_mode::none ) {
        std::cout << "none" << std::endl; }
    else if( m_URR_mode == Transporting::URR_mode::pdfs ) {
        std::cout << "pdfs" << std::endl; }
    else if( m_URR_mode == Transporting::URR_mode::ACE_URR_probabilityTables ) {
        std::cout << "ACE style protability tables" << std::endl; }
    else {
        std::cout << "Oops, need to code into print method." << std::endl;
    }

    std::cout << a_indent << "# Index of reactions in URR region:";
    for( auto reactionIter = m_reactionsInURR_region.begin( ); reactionIter != m_reactionsInURR_region.end( ); ++reactionIter ) {
        std::cout << indent2 << LUPI::Misc::argumentsToString( a_iFormat.c_str( ), *reactionIter ) << std::endl;
    }
    std::cout << std::endl;

    int reactionIndex = 0;
    std::cout << std::endl;
    std::cout << a_indent << "# Number of reactions = " << m_reactionCrossSections.size( ) << std::endl;
    for( auto iter = m_reactionCrossSections.begin( ); iter != m_reactionCrossSections.end( ); ++iter, ++reactionIndex ) {
        Reaction const *reaction = a_protareSingle->reaction( reactionIndex );

        std::cout << a_indent << "# Reaction number " << reactionIndex << std::endl;
        std::cout << a_indent << "# Reaction label: " << reaction->label( ).c_str( ) << std::endl;
        (*iter)->print( a_protareSingle, a_indent, a_iFormat, a_energyFormat, a_dFormat );
    }
}

/*
============================================================
=========== HeatedCrossSectionsContinuousEnergy ============
============================================================
*/
LUPI_HOST_DEVICE HeatedCrossSectionsContinuousEnergy::HeatedCrossSectionsContinuousEnergy( ) :
        m_temperatures( ),
        m_thresholds( ),
        m_heatedCrossSections( ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE HeatedCrossSectionsContinuousEnergy::~HeatedCrossSectionsContinuousEnergy( ) {

    for( Vector<HeatedCrossSectionContinuousEnergy *>::const_iterator iter = m_heatedCrossSections.begin( ); iter != m_heatedCrossSections.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Fills in *this* with the requested temperature data.
 *
 * @param a_smr                         [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_setupInfo                   [in]    Used internally when constructing a Protare to pass information to other components.
 * @param a_settings                    [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles                   [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 * @param a_domainHash                  [in]    The hash data used when looking up a cross section.
 * @param a_temperatureInfos            [in]    The list of temperatures to use.
 * @param a_reactions                   [in]    The list of reactions to use.
 * @param a_orphanProducts              [in]    The list of orphan products to use.
 * @param a_fixedGrid                   [in]    If true, the specified fixed grid is used; otherwise, grid in the file is used.
 * @param a_zeroReactions               [in]    Special case where no reaction in a protare is wanted so the first one is used but its cross section is set to 0.0 at all energies.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionsContinuousEnergy::update( LUPI_maybeUnused LUPI::StatusMessageReporting &a_smr, SetupInfo &a_setupInfo, 
                Transporting::MC const &a_settings, GIDI::Transporting::Particles const &a_particles, DomainHash const &a_domainHash, 
                GIDI::Styles::TemperatureInfos const &a_temperatureInfos, std::vector<GIDI::Reaction const *> const &a_reactions, 
                std::vector<GIDI::Reaction const *> const &a_orphanProducts, bool a_fixedGrid, bool a_zeroReactions ) {

    m_temperatures.reserve( a_temperatureInfos.size( ) );
    m_heatedCrossSections.reserve( a_temperatureInfos.size( ) );

    for( GIDI::Styles::TemperatureInfos::const_iterator iter = a_temperatureInfos.begin( ); iter != a_temperatureInfos.end( ); ++iter ) {
        m_temperatures.push_back( iter->temperature( ).value( ) );
        m_heatedCrossSections.push_back( new HeatedCrossSectionContinuousEnergy( a_setupInfo, a_settings, a_particles, a_domainHash, *iter, 
                a_reactions, a_orphanProducts, a_fixedGrid, a_zeroReactions ) );
    }

    m_thresholds.resize( m_heatedCrossSections[0]->numberOfReactions( ) );
    for( int i1 = 0; i1 < m_heatedCrossSections[0]->numberOfReactions( ); ++i1 ) m_thresholds[i1] = m_heatedCrossSections[0]->threshold( i1 );
}

/* *********************************************************************************************************//**
 * Returns the cross section for target temperature *a_temperature* and projectile energy *a_energy*.
 *
 * @param a_URR_protareInfos    [in]    URR information.
 * @param a_URR_index           [in]    If not negative, specifies the index in *a_URR_protareInfos*.
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_sampling            [in]    If *true* the cross section is to be used for sampling, otherwise, just for looking up.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsContinuousEnergy::crossSection( URR_protareInfos const &a_URR_protareInfos, int a_URR_index, 
                int a_hashIndex, double a_temperature, double a_energy, bool a_sampling ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double cross_section;

    if( a_temperature <= m_temperatures[0] ) { 
        cross_section = m_heatedCrossSections[0]->crossSection( a_URR_protareInfos, a_URR_index, a_hashIndex, a_energy, a_sampling ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        cross_section = m_heatedCrossSections.back( )->crossSection( a_URR_protareInfos, a_URR_index, a_hashIndex, a_energy, a_sampling ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        cross_section = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->crossSection( a_URR_protareInfos, a_URR_index, a_hashIndex, a_energy, a_sampling )
                            + fraction * m_heatedCrossSections[i1]->crossSection( a_URR_protareInfos, a_URR_index, a_hashIndex, a_energy, a_sampling );
    }

    return( cross_section );
}

/* *********************************************************************************************************//**
 * Adds the energy dependent, total cross section corresponding to the temperature *a_temperature* multiplied by *a_userFactor* to *a_crossSectionVector*.
 * This function only works for fixed-grid data.
 *
 * @param   a_temperature               [in]        Specifies the temperature of the material.
 * @param   a_userFactor                [in]        User factor which all cross sections are multiplied by.
 * @param   a_numberAllocated           [in]        The length of memory allocated for *a_crossSectionVector*.
 * @param   a_crossSectionVector        [in/out]    The energy dependent, total cross section to add cross section data to.
 ***********************************************************************************************************/
 
LUPI_HOST_DEVICE void HeatedCrossSectionsContinuousEnergy::crossSectionVector( double a_temperature, double a_userFactor, 
                std::size_t a_numberAllocated, double *a_crossSectionVector ) const {

    int number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    int index1 = 0, index2 = 0;
    double fraction = 0.0;

    if( a_temperature <= m_temperatures[0] ) { 
        }
    else if( a_temperature >= m_temperatures.back( ) ) {
        index1 = index2 = number_of_temperatures - 1;
        fraction = 1.0; }
    else {
        for( ; index2 < number_of_temperatures; ++index2 ) if( a_temperature < m_temperatures[index2] ) break;
        index1 = index2 - 1;
        fraction = ( a_temperature - m_temperatures[index1] ) / ( m_temperatures[index2] - m_temperatures[index1] );
    }

    Vector<MCGIDI_FLOAT> &totalCrossSection1 = m_heatedCrossSections[index1]->totalCrossSection( );
    Vector<MCGIDI_FLOAT> &totalCrossSection2 = m_heatedCrossSections[index2]->totalCrossSection( );
    std::size_t size = totalCrossSection1.size( );
    double factor1 = a_userFactor * ( 1.0 - fraction ), factor2 = a_userFactor * fraction;

    if( a_numberAllocated < totalCrossSection1.size( ) ) LUPI_THROW( "HeatedCrossSectionsContinuousEnergy::crossSectionVector: a_numberAllocated too small." );
    for( std::size_t i1 = 0; i1 < size; ++i1 ) {
        a_crossSectionVector[i1] += factor1 * totalCrossSection1[i1] + factor2 * totalCrossSection2[i1];
    }
}

/* *********************************************************************************************************//**
 * Returns the total cross section as a GIDI::Functions::XYs1d instance.
 *
 * @param a_temperature         [in]    Specifies the temperature requested for the cross section.
 * 
 * @returns                             A GIDI::Functions::XYs1d instance.
 ***********************************************************************************************************/

LUPI_HOST GIDI::Functions::XYs1d HeatedCrossSectionsContinuousEnergy::crossSectionAsGIDI_XYs1d( double a_temperature ) const {

    GIDI::Functions::XYs1d crossSection1;

    if( a_temperature <= m_temperatures[0] ) {
        crossSection1 = m_heatedCrossSections[0]->crossSectionAsGIDI_XYs1d( ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        crossSection1 = m_heatedCrossSections.back( )->crossSectionAsGIDI_XYs1d( ); }
    else {
        int number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
        int i1 = 0;
        for( ; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        crossSection1 = m_heatedCrossSections[i1-1]->crossSectionAsGIDI_XYs1d( );
        crossSection1 *= ( 1. - fraction );
        GIDI::Functions::XYs1d crossSection2( m_heatedCrossSections[i1-1]->crossSectionAsGIDI_XYs1d( ) );
        crossSection2 *= fraction;
        crossSection1 += crossSection2;
    }

    return( crossSection1 );
}

/* *********************************************************************************************************//**
 * Returns the reaction's cross section as GIDI::Functions::XYs1d instance.
 *
 * @param a_reactionIndex       [in]    Specifies the indexs of the reaction whose cross section is requested.
 * @param a_temperature         [in]    Specifies the temperature requested for the cross section.
 *
 * @returns                             A GIDI::Functions::XYs1d instance.
 ***********************************************************************************************************/

LUPI_HOST GIDI::Functions::XYs1d HeatedCrossSectionsContinuousEnergy::reactionCrossSectionAsGIDI_XYs1d( int a_reactionIndex, 
                double a_temperature ) const {

    GIDI::Functions::XYs1d crossSection1;

    if( a_temperature <= m_temperatures[0] ) {
        crossSection1 = m_heatedCrossSections[0]->reactionCrossSectionAsGIDI_XYs1d( a_reactionIndex ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        crossSection1 = m_heatedCrossSections.back( )->reactionCrossSectionAsGIDI_XYs1d( a_reactionIndex ); }
    else {
        int number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
        int i1 = 0;
        for( ; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        crossSection1 = m_heatedCrossSections[i1-1]->reactionCrossSectionAsGIDI_XYs1d( a_reactionIndex );
        crossSection1 *= ( 1. - fraction );
        GIDI::Functions::XYs1d crossSection2( m_heatedCrossSections[i1-1]->reactionCrossSectionAsGIDI_XYs1d( a_reactionIndex ) );
        crossSection2 *= fraction;
        crossSection1 += crossSection2;
    }

    return( crossSection1 );
}

/*
=========================================================
*/
LUPI_HOST_DEVICE double HeatedCrossSectionsContinuousEnergy::reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_URR_index, 
                int a_hashIndex, double a_temperature, double a_energy, bool a_sampling ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double cross_section;

    if( a_temperature <= m_temperatures[0] ) {
        cross_section = m_heatedCrossSections[0]->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_URR_index, a_hashIndex, a_energy, a_sampling ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        cross_section = m_heatedCrossSections.back( )->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_URR_index, a_hashIndex, a_energy, a_sampling ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        cross_section = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_URR_index, a_hashIndex, a_energy, a_sampling )
                            + fraction * m_heatedCrossSections[i1]->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_URR_index, a_hashIndex, a_energy, a_sampling );
    }

    return( cross_section );
}

/*
=========================================================
*/
LUPI_HOST_DEVICE double HeatedCrossSectionsContinuousEnergy::reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_URR_index,
                double a_temperature, double a_energy_in ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double cross_section;

    if( a_temperature <= m_temperatures[0] ) {
        cross_section = m_heatedCrossSections[0]->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_URR_index, a_energy_in ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        cross_section = m_heatedCrossSections.back( )->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_URR_index, a_energy_in ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        cross_section = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_URR_index, a_energy_in )
                            + fraction * m_heatedCrossSections[i1]->reactionCrossSection( a_reactionIndex, a_URR_protareInfos, a_URR_index, a_energy_in );
    }

    return( cross_section );
}

/* *********************************************************************************************************//**
 * Returns the deposition energy for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsContinuousEnergy::depositionEnergy( int a_hashIndex, double a_temperature, double a_energy ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double deposition_energy;

    if( a_temperature <= m_temperatures[0] ) { 
        deposition_energy = m_heatedCrossSections[0]->depositionEnergy( a_hashIndex, a_energy ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        deposition_energy = m_heatedCrossSections.back( )->depositionEnergy( a_hashIndex, a_energy ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        deposition_energy = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->depositionEnergy( a_hashIndex, a_energy )
                            + fraction * m_heatedCrossSections[i1]->depositionEnergy( a_hashIndex, a_energy );
    }

    return( deposition_energy );
}

/* *********************************************************************************************************//**
 * Returns the deposition momentum for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsContinuousEnergy::depositionMomentum( int a_hashIndex, double a_temperature, double a_energy ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double deposition_momentum;

    if( a_temperature <= m_temperatures[0] ) {
        deposition_momentum = m_heatedCrossSections[0]->depositionMomentum( a_hashIndex, a_energy ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        deposition_momentum = m_heatedCrossSections.back( )->depositionMomentum( a_hashIndex, a_energy ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        deposition_momentum = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->depositionMomentum( a_hashIndex, a_energy )
                            + fraction * m_heatedCrossSections[i1]->depositionMomentum( a_hashIndex, a_energy );
    }

    return( deposition_momentum );
}

/* *********************************************************************************************************//**
 * Returns the production momentum for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsContinuousEnergy::productionEnergy( int a_hashIndex, double a_temperature, double a_energy ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double production_energy;

    if( a_temperature <= m_temperatures[0] ) {
        production_energy = m_heatedCrossSections[0]->productionEnergy( a_hashIndex, a_energy ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        production_energy = m_heatedCrossSections.back( )->productionEnergy( a_hashIndex, a_energy ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        production_energy = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->productionEnergy( a_hashIndex, a_energy )
                            + fraction * m_heatedCrossSections[i1]->productionEnergy( a_hashIndex, a_energy );
    }

    return( production_energy );
}

/* *********************************************************************************************************//**
 * Returns the gain for particle with index *a_particleIndex* for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_particleIndex       [in]    The index of the particle whose gain is requested.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsContinuousEnergy::gain( int a_hashIndex, double a_temperature, double a_energy, int a_particleIndex ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double production_energy;

    if( a_temperature <= m_temperatures[0] ) {
        production_energy = m_heatedCrossSections[0]->gain( a_hashIndex, a_energy, a_particleIndex ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        production_energy = m_heatedCrossSections.back( )->gain( a_hashIndex, a_energy, a_particleIndex ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        production_energy = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->gain( a_hashIndex, a_energy, a_particleIndex )
                            + fraction * m_heatedCrossSections[i1]->gain( a_hashIndex, a_energy, a_particleIndex );
    }

    return( production_energy );
}

/* *********************************************************************************************************//**
 * Returns the gain for particle with intid *a_particleIntid* for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_hashIndex           [in]    Specifies projectile energy hash index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_particleIntid       [in]    The intid of the particle whose gain is requested.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsContinuousEnergy::gainViaIntid( int a_hashIndex, double a_temperature, double a_energy, int a_particleIntid ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double production_energy;

    if( a_temperature <= m_temperatures[0] ) {
        production_energy = m_heatedCrossSections[0]->gainViaIntid( a_hashIndex, a_energy, a_particleIntid ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        production_energy = m_heatedCrossSections.back( )->gainViaIntid( a_hashIndex, a_energy, a_particleIntid ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        production_energy = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->gainViaIntid( a_hashIndex, a_energy, a_particleIntid )
                            + fraction * m_heatedCrossSections[i1]->gainViaIntid( a_hashIndex, a_energy, a_particleIntid );
    }

    return( production_energy );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *
 * @param a_particleIndex       [in]    The index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionsContinuousEnergy::setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {

    for( auto iter = m_heatedCrossSections.begin( ); iter != m_heatedCrossSections.end( ); ++iter ) (*iter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *
 * @param a_particleIntid       [in]    The intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionsContinuousEnergy::setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {

    for( auto iter = m_heatedCrossSections.begin( ); iter != m_heatedCrossSections.end( ); ++iter ) (*iter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void HeatedCrossSectionsContinuousEnergy::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_VECTOR_DOUBLE( m_temperatures, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_thresholds, a_buffer, a_mode );

    std::size_t vectorSize = m_heatedCrossSections.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_heatedCrossSections.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_heatedCrossSections.internalSize();
    for( std::size_t memberIndex = 0; memberIndex < vectorSize; ++memberIndex ) {
        if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
            if( a_buffer.m_placement != nullptr ) {
                m_heatedCrossSections[memberIndex] = new(a_buffer.m_placement) HeatedCrossSectionContinuousEnergy;
                a_buffer.incrementPlacement( sizeof( HeatedCrossSectionContinuousEnergy ) ); }
            else {
                m_heatedCrossSections[memberIndex] = new HeatedCrossSectionContinuousEnergy;
            }
        }
        if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof( HeatedCrossSectionContinuousEnergy ) );
        }
        m_heatedCrossSections[memberIndex]->serialize( a_buffer, a_mode );
    }
}

/* *********************************************************************************************************//**
 * Print to *std::cout* the content of *this*. This is mainly meant for debugging.
 *
 * @param a_protareSingle       [in]    The GIDI::ProtareSingle instance that contains *this*.
 * @param a_indent              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_iFormat             [in]    C printf format specifier for any interger that is printed (e.g., "%3d").
 * @param a_energyFormat        [in]    C printf format specifier for any interger that is printed (e.g., "%20.12e").
 * @param a_dFormat             [in]    C printf format specifier for any interger that is printed (e.g., "%14.7e").
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionsContinuousEnergy::print( ProtareSingle const *a_protareSingle, std::string const &a_indent,
                std::string const &a_iFormat, std::string const &a_energyFormat, std::string const &a_dFormat ) const {

    std::cout << "Temperatures present:" << std::endl;
    for( auto temperatureIter = m_temperatures.begin( ); temperatureIter != m_temperatures.end( ); ++temperatureIter ) {
        std::cout << LUPI::Misc::argumentsToString( a_dFormat.c_str( ), *temperatureIter ) << std::endl;
    }
    std::cout << "Reaction thresholds:" << std::endl;
    for( auto reactionIter = m_thresholds.begin( ); reactionIter != m_thresholds.end( ); ++reactionIter ) {
        std::cout << LUPI::Misc::argumentsToString( a_dFormat.c_str( ), *reactionIter ) << std::endl;
    }
    for( auto heatedCrossSections = m_heatedCrossSections.begin( ); heatedCrossSections != m_heatedCrossSections.end( ); ++heatedCrossSections ) {
        (*heatedCrossSections)->print( a_protareSingle, a_indent, a_iFormat, a_energyFormat, a_dFormat );
    }
}

/*! \class MultiGroupGain
 * This class store a particles index and gain for a protare.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE MultiGroupGain::MultiGroupGain( ) :
        m_particleIntid( -1 ),
        m_particleIndex( -1 ),
        m_userParticleIndex( -1 ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST MultiGroupGain::MultiGroupGain( int a_particleIntid, int a_particleIndex, GIDI::Vector const &a_gain ) :
        m_particleIntid( a_particleIntid ),
        m_particleIndex( a_particleIndex ),
        m_userParticleIndex( -1 ),
        m_gain( GIDI_VectorDoublesToMCGIDI_VectorDoubles( a_gain ) ) {

}

/* *********************************************************************************************************//**
 * @param a_multiGroupGain      [in]    The **MultiGroupGain** whose contents are to be copied.
 ***********************************************************************************************************/

LUPI_HOST MultiGroupGain &MultiGroupGain::operator=( MultiGroupGain const &a_multiGroupGain ) {

    m_particleIntid = a_multiGroupGain.particleIntid( );
    m_particleIndex = a_multiGroupGain.particleIndex( );
    m_userParticleIndex = a_multiGroupGain.userParticleIndex( );
    m_gain = a_multiGroupGain.gain( );

    return( *this );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void MultiGroupGain::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_INT( m_particleIntid, a_buffer, a_mode );
    DATA_MEMBER_INT( m_particleIndex, a_buffer, a_mode );
    DATA_MEMBER_INT( m_userParticleIndex, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_gain, a_buffer, a_mode );
}

/* *********************************************************************************************************//**
 * This method writes the multi-group data.
 *
 * @param a_file                [in]    The buffer to read or write data to depending on *a_mode*.
 ***********************************************************************************************************/

LUPI_HOST void MultiGroupGain::write( FILE *a_file ) const {

    std::string buffer = LUPI::Misc::argumentsToString( "Gain for particle %d (%d, %d)", m_particleIntid, m_particleIndex, m_userParticleIndex );
    writeVector( a_file, buffer, 0, m_gain );
}

/*
============================================================
=========== HeatedReactionCrossSectionMultiGroup ===========
============================================================
*/
LUPI_HOST_DEVICE HeatedReactionCrossSectionMultiGroup::HeatedReactionCrossSectionMultiGroup( ) {

}
/*
=========================================================
*/
LUPI_HOST HeatedReactionCrossSectionMultiGroup::HeatedReactionCrossSectionMultiGroup( SetupInfo &a_setupInfo, LUPI_maybeUnused Transporting::MC const &a_settings, 
                int a_offset, std::vector<double> const &a_crossSection, double a_threshold ) :
        m_threshold( a_threshold ),
        m_offset( a_offset ),
        m_crossSections( a_crossSection ),
        m_augmentedThresholdCrossSection( 0.0 ) {

    Vector<double> const &boundaries = a_setupInfo.m_protare.projectileMultiGroupBoundariesCollapsed( );

    if( ( a_offset > 0 ) && ( boundaries[a_offset] < a_threshold ) ) {  // This uses the linear rejection above threshold in the group m_offset.
        if( ( boundaries[a_offset] < a_threshold ) && ( a_threshold < boundaries[a_offset+1] ) ) {
            m_augmentedThresholdCrossSection = m_crossSections[0] * ( boundaries[a_offset+1] + a_threshold - 2.0 * boundaries[a_offset] ) / ( boundaries[a_offset+1] - a_threshold ); }
        else {
            m_crossSections[0] = 0.0;
        }
    }
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void HeatedReactionCrossSectionMultiGroup::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_DOUBLE( m_threshold, a_buffer, a_mode  );
    DATA_MEMBER_INT( m_offset, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_crossSections, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_augmentedThresholdCrossSection, a_buffer, a_mode  );
}

/* *********************************************************************************************************//**
 * This method writes the multi-group data.
 *
 * @param a_file                [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_reactionIndex       [in]    The index of the reaction.
 ***********************************************************************************************************/

LUPI_HOST void HeatedReactionCrossSectionMultiGroup::write( FILE *a_file, int a_reactionIndex ) const {

    std::string buffer = LUPI::Misc::argumentsToString( "Reaction cross section (%3d)", a_reactionIndex );
    writeVector( a_file, buffer, m_offset, m_crossSections );
}

/*
============================================================
============= HeatedCrossSectionMultiGroup ================
============================================================
*/
LUPI_HOST_DEVICE HeatedCrossSectionMultiGroup::HeatedCrossSectionMultiGroup( ) {

}
/*
=========================================================
*/
LUPI_HOST HeatedCrossSectionMultiGroup::HeatedCrossSectionMultiGroup( LUPI::StatusMessageReporting &a_smr, GIDI::ProtareSingle const &a_protare, 
                SetupInfo &a_setupInfo, Transporting::MC const &a_settings, GIDI::Styles::TemperatureInfo const &a_temperatureInfo, 
                GIDI::Transporting::Particles const &a_particles, std::vector<GIDI::Reaction const *> const &a_reactions, LUPI_maybeUnused std::string const &a_label,
                bool a_zeroReactions, GIDI::ExcludeReactionsSet const &a_reactionsToExclude ) :
        m_totalCrossSection( ),
        m_augmentedCrossSection( ),
        m_reactionCrossSections( ) {

    GIDI::Transporting::Mode transportMode = GIDI::Transporting::Mode::multiGroup;
    if( ( a_settings.upscatterModel( ) == Sampling::Upscatter::Model::B ) || ( a_settings.upscatterModel( ) == Sampling::Upscatter::Model::BSnLimits ) 
            || ( a_settings.upscatterModel( ) == Sampling::Upscatter::Model::A ) )
        transportMode = GIDI::Transporting::Mode::multiGroupWithSnElasticUpScatter;
    GIDI::Transporting::MG multi_group_settings( a_settings.projectileID( ), transportMode, a_settings.delayedNeutrons( ) );
    multi_group_settings.setThrowOnError( a_settings.throwOnError( ) );

    GIDI::Axes axes;
    std::vector<double> dummy;
    GIDI::Vector totalCrossSection;
    GIDI::Vector vector;

    m_reactionCrossSections.reserve( a_reactions.size( ) );
    int index = 0;                                      // Only used for debugging.
    GIDI::Transporting::MG MG_settings( a_settings.projectileID( ), GIDI::Transporting::Mode::multiGroup, a_settings.delayedNeutrons( ) );
    MG_settings.setThrowOnError( a_settings.throwOnError( ) );

    for( std::vector<GIDI::Reaction const *>::const_iterator reactionIter = a_reactions.begin( ); reactionIter != a_reactions.end( ); ++reactionIter, ++index ) {
        GIDI::Vector crossSectionVector = (*reactionIter)->multiGroupCrossSection( a_smr, MG_settings, a_temperatureInfo );

        vector = GIDI::collapse( crossSectionVector, a_settings, a_particles, 0.0 );

        std::size_t start = 0;
        for( ; start < vector.size( ); ++start ) {
            if( vector[start] != 0.0 ) break;
        }
        checkZeroReaction( vector, a_zeroReactions );
        std::vector<double> data;
        for( std::size_t i1 = start; i1 < vector.size( ); ++i1 ) data.push_back( vector[i1] );
        int offset = static_cast<int>( start );
    
        m_reactionCrossSections.push_back( new HeatedReactionCrossSectionMultiGroup( a_setupInfo, a_settings, offset, data, 
                (*reactionIter)->crossSectionThreshold( ) ) );

        totalCrossSection += vector;
    }

    m_totalCrossSection = totalCrossSection.data( );

    m_augmentedCrossSection.resize( totalCrossSection.size( ) );
    for( std::size_t i1 = 0; i1 < m_augmentedCrossSection.size( ); ++i1 ) m_augmentedCrossSection[i1] = 0;
    for( std::size_t i1 = 0; i1 < m_reactionCrossSections.size( ); ++i1 )
        m_augmentedCrossSection[m_reactionCrossSections[i1]->offset( )] += m_reactionCrossSections[i1]->augmentedThresholdCrossSection( );


    vector = a_protare.multiGroupDepositionEnergy( a_smr, multi_group_settings, a_temperatureInfo, a_particles, a_reactionsToExclude );
    vector = collapseAndcheckZeroReaction( vector, a_settings, a_particles, 0.0, a_zeroReactions );
    m_depositionEnergy = GIDI_VectorDoublesToMCGIDI_VectorDoubles( vector );

    vector = a_protare.multiGroupDepositionMomentum( a_smr, multi_group_settings, a_temperatureInfo, a_particles, a_reactionsToExclude );
    vector = collapseAndcheckZeroReaction( vector, a_settings, a_particles, 0.0, a_zeroReactions );
    m_depositionMomentum = GIDI_VectorDoublesToMCGIDI_VectorDoubles( vector );

    vector = a_protare.multiGroupQ( a_smr, multi_group_settings, a_temperatureInfo, true, true, a_reactionsToExclude );
    vector = collapseAndcheckZeroReaction( vector, a_settings, a_particles, 0.0, a_zeroReactions );
    m_productionEnergy = GIDI_VectorDoublesToMCGIDI_VectorDoubles( vector );

    std::map<std::string, GIDI::Transporting::Particle> particles = a_particles.particles( );
    m_gains.resize( particles.size( ) );
    int i1 = 0;
    for( std::map<std::string, GIDI::Transporting::Particle>::const_iterator particle = particles.begin( ); particle != particles.end( ); ++particle, ++i1 ) {
        int particleIntid = a_setupInfo.m_particleIntids[particle->first];
        int particleIndex = a_setupInfo.m_particleIndices[particle->first];

        vector = a_protare.multiGroupGain( a_smr, multi_group_settings, a_temperatureInfo, particle->first, a_reactionsToExclude );
        vector = collapseAndcheckZeroReaction( vector, a_settings, a_particles, 0.0, a_zeroReactions );
        m_gains[i1] = new MultiGroupGain( particleIntid, particleIndex, vector );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE HeatedCrossSectionMultiGroup::~HeatedCrossSectionMultiGroup( ) {

    for( auto iter = m_reactionCrossSections.begin( ); iter < m_reactionCrossSections.end( ); ++iter ) delete *iter;
    for( auto iter = m_gains.begin( ); iter < m_gains.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Returns the multi-group cross section.
 *
 * @param a_hashIndex           [in]    The multi-group index.
 * @param a_sampling            [in]    Fix me.
 *
 * @return                              A vector of the length of the number of multi-group groups.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionMultiGroup::crossSection( int a_hashIndex, bool a_sampling ) const { 

    double crossSection2 = m_totalCrossSection[a_hashIndex];

    if( a_sampling ) crossSection2 += m_augmentedCrossSection[a_hashIndex];

    return( crossSection2 );
}

/* *********************************************************************************************************//**
 * Returns the multi-group gain for particle with index *a_particleIndex*. If no particle is found, a Vector of all 0's is returned.
 *
 * @param a_particleIndex       [in]    The index of the particle whose gain is to be returned.
 * @param a_hashIndex           [in]    The multi-group index.
 *
 * @return                              A vector of the length of the number of multi-group groups.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionMultiGroup::gain( int a_particleIndex, int a_hashIndex ) const {

    for( std::size_t i1 = 0; i1 < m_gains.size( ); ++i1 ) {
        if( a_particleIndex == m_gains[i1]->particleIndex( ) ) return( m_gains[i1]->gain( a_hashIndex ) );
    }

    return( 0.0 );
}

/* *********************************************************************************************************//**
 * Returns the multi-group gain for particle with intid *a_particleIntid*. If no particle is found, a Vector of all 0's is returned.
 *
 * @param a_particleIntid       [in]    The intid of the particle whose gain is to be returned.
 * @param a_hashIndex           [in]    The multi-group index.
 *
 * @return                              A vector of the length of the number of multi-group groups.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionMultiGroup::gainViaIntid( int a_particleIntid, int a_hashIndex ) const {

    for( std::size_t i1 = 0; i1 < m_gains.size( ); ++i1 ) {
        if( a_particleIntid == m_gains[i1]->particleIntid( ) ) return( m_gains[i1]->gain( a_hashIndex ) );
    }

    return( 0.0 );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIntid*.
 *
 * @param a_particleIndex       [in]    The index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionMultiGroup::setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {

    for( auto iter = m_gains.begin( ); iter != m_gains.end( ); ++iter ) (*iter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *
 * @param a_particleIntid       [in]    The intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionMultiGroup::setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {

    for( auto iter = m_gains.begin( ); iter != m_gains.end( ); ++iter ) (*iter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void HeatedCrossSectionMultiGroup::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_VECTOR_DOUBLE( m_totalCrossSection, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_augmentedCrossSection, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_depositionEnergy, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_depositionMomentum, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_productionEnergy, a_buffer, a_mode );

    std::size_t vectorSize = m_reactionCrossSections.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_reactionCrossSections.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_reactionCrossSections.internalSize();
    for( std::size_t memberIndex = 0; memberIndex < vectorSize; ++memberIndex ) {
        if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
            if( a_buffer.m_placement != nullptr ) {
                m_reactionCrossSections[memberIndex] = new(a_buffer.m_placement) HeatedReactionCrossSectionMultiGroup;
                a_buffer.incrementPlacement( sizeof( HeatedReactionCrossSectionMultiGroup ) ); }
            else {
                m_reactionCrossSections[memberIndex] = new HeatedReactionCrossSectionMultiGroup;
            }
        }
        if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof( HeatedReactionCrossSectionMultiGroup ) );
        }
        m_reactionCrossSections[memberIndex]->serialize( a_buffer, a_mode );
    }

    vectorSize = m_gains.size( );
    vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_gains.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_gains.internalSize();
    for( std::size_t memberIndex = 0; memberIndex < vectorSize; ++memberIndex ) {
        if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
            if( a_buffer.m_placement != nullptr ) {
                m_gains[memberIndex] = new(a_buffer.m_placement) MultiGroupGain;
                a_buffer.incrementPlacement( sizeof( MultiGroupGain ) ); }
            else {
                m_gains[memberIndex] = new MultiGroupGain;
            }
        }
        if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof( MultiGroupGain ) );
        }
        m_gains[memberIndex]->serialize( a_buffer, a_mode );
    }
}

/* *********************************************************************************************************//**
 * This method writes the multi-group data.
 *
 * @param a_file                [in]    The buffer to read or write data to depending on *a_mode*.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionMultiGroup::write( FILE *a_file ) const {

    writeVector( a_file, "Total cross section", 0, m_totalCrossSection );
    writeVector( a_file, "Augmented cross section", 0, m_augmentedCrossSection );
    writeVector( a_file, "Deposition energy", 0, m_depositionEnergy );
    writeVector( a_file, "Deposition momentum", 0, m_depositionMomentum );
    writeVector( a_file, "Production energy", 0, m_productionEnergy );

    for( auto iter = m_gains.begin( ); iter != m_gains.end( ); ++iter ) (*iter)->write( a_file );
    int reactionIndex = 0;
    for( Vector<HeatedReactionCrossSectionMultiGroup *>::const_iterator iter = m_reactionCrossSections.begin( ); iter < m_reactionCrossSections.end( ); ++iter ) {
        (*iter)->write( a_file, reactionIndex );
        ++reactionIndex;
    }
}

/*! \class Protare
 * Base class for the protare sub-classes.
 */

/* *********************************************************************************************************//**
 * Generic constructor.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE HeatedCrossSectionsMultiGroup::HeatedCrossSectionsMultiGroup( ) {

}

/* *********************************************************************************************************//**
 * Generic constructor.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE HeatedCrossSectionsMultiGroup::~HeatedCrossSectionsMultiGroup( ) {

    for( Vector<HeatedCrossSectionMultiGroup *>::const_iterator iter = m_heatedCrossSections.begin( ); iter != m_heatedCrossSections.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Fills in *this* with the requested temperature data.
 *
 * @param a_smr                         [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_protare                     [in]    The GIDI::Protare used to constuct the Protare that *this* is a part of.
 * @param a_setupInfo                   [in]    Used internally when constructing a Protare to pass information to other components.
 * @param a_settings                    [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles                   [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 * @param a_temperatureInfos            [in]    The list of temperatures to use.
 * @param a_reactions                   [in]    The list of reactions to use.
 * @param a_orphanProducts              [in]    The list of orphan products to use.
 * @param a_zeroReactions               [in]    Special case where no reaction in a protare is wanted so the first one is used but its cross section is set to 0.0 at all energies.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionsMultiGroup::update( LUPI::StatusMessageReporting &a_smr, GIDI::ProtareSingle const &a_protare, 
                SetupInfo &a_setupInfo, Transporting::MC const &a_settings, GIDI::Transporting::Particles const &a_particles, 
                GIDI::Styles::TemperatureInfos const &a_temperatureInfos, std::vector<GIDI::Reaction const *> const &a_reactions, 
                LUPI_maybeUnused std::vector<GIDI::Reaction const *> const &a_orphanProducts, bool a_zeroReactions,
                GIDI::ExcludeReactionsSet const &a_reactionsToExclude ) {

    m_temperatures.reserve( a_temperatureInfos.size( ) );
    m_heatedCrossSections.reserve( a_temperatureInfos.size( ) );

    for( GIDI::Styles::TemperatureInfos::const_iterator iter = a_temperatureInfos.begin( ); iter != a_temperatureInfos.end( ); ++iter ) {
        m_temperatures.push_back( iter->temperature( ).value( ) );
        m_heatedCrossSections.push_back( new HeatedCrossSectionMultiGroup( a_smr, a_protare, a_setupInfo, a_settings, *iter, a_particles, 
                a_reactions, iter->heatedMultiGroup( ), a_zeroReactions, a_reactionsToExclude ) );
    }

    m_thresholds.resize( m_heatedCrossSections[0]->numberOfReactions( ) );
    for( int i1 = 0; i1 < m_heatedCrossSections[0]->numberOfReactions( ); ++i1 ) m_thresholds[i1] = m_heatedCrossSections[0]->threshold( i1 );

    m_multiGroupThresholdIndex.resize( m_heatedCrossSections[0]->numberOfReactions( ) );
    for( int i1 = 0; i1 < m_heatedCrossSections[0]->numberOfReactions( ); ++i1 ) {
        m_multiGroupThresholdIndex[i1] = -1;
        if( m_thresholds[i1] > 0 ) m_multiGroupThresholdIndex[i1] = m_heatedCrossSections[0]->thresholdOffset( i1 );
    }

    m_projectileMultiGroupBoundariesCollapsed = a_setupInfo.m_protare.projectileMultiGroupBoundariesCollapsed( );
}

/* *********************************************************************************************************//**
 * Returns the total multi-group cross section for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_hashIndex           [in]    The multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_sampling            [in]    Used for multi-group look up. If *true*, use augmented cross sections.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsMultiGroup::crossSection( int a_hashIndex, double a_temperature, bool a_sampling ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double cross_section;

    if( a_temperature <= m_temperatures[0] ) {
        cross_section = m_heatedCrossSections[0]->crossSection( a_hashIndex, a_sampling ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        cross_section = m_heatedCrossSections.back( )->crossSection( a_hashIndex, a_sampling ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );

        cross_section = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->crossSection( a_hashIndex, a_sampling )
                        + fraction * m_heatedCrossSections[i1]->crossSection( a_hashIndex, a_sampling );
    }

    return( cross_section );
}

/* *********************************************************************************************************//**
 * Adds the energy dependent, total cross section corresponding to the temperature *a_temperature* multiplied by *a_userFactor* to *a_crossSectionVector*.
 *
 * @param   a_temperature                   [in]        Specifies the temperature of the material.
 * @param   a_userFactor                    [in]        User factor which all cross sections are multiplied by.
 * @param   a_numberAllocated               [in]        The length of memory allocated for *a_crossSectionVector*.
 * @param   a_crossSectionVector            [in/out]    The energy dependent, total cross section to add cross section data to.
 ***********************************************************************************************************/
 
LUPI_HOST_DEVICE void HeatedCrossSectionsMultiGroup::crossSectionVector( double a_temperature, double a_userFactor, 
                std::size_t a_numberAllocated, double *a_crossSectionVector ) const {

    int number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    int index1 = 0, index2 = 0;
    double fraction = 0.0;

    if( a_temperature <= m_temperatures[0] ) { 
        }
    else if( a_temperature >= m_temperatures.back( ) ) {
        index1 = index2 = number_of_temperatures - 1;
        fraction = 1.0; }
    else {
        for( ; index2 < number_of_temperatures; ++index2 ) if( a_temperature < m_temperatures[index2] ) break;
        index1 = index2 - 1;
        fraction = ( a_temperature - m_temperatures[index1] ) / ( m_temperatures[index2] - m_temperatures[index1] );
    }

    Vector<double> &totalCrossSection1 = m_heatedCrossSections[index1]->totalCrossSection( );
    Vector<double> &totalCrossSection2 = m_heatedCrossSections[index2]->totalCrossSection( );
    std::size_t size = totalCrossSection1.size( );
    double factor1 = a_userFactor * ( 1.0 - fraction ), factor2 = a_userFactor * fraction;

    if( a_numberAllocated < totalCrossSection1.size( ) ) LUPI_THROW( "HeatedCrossSectionsMultiGroup::crossSectionVector: a_numberAllocated too small." );
    for( std::size_t i1 = 0; i1 < size; ++i1 ) {
        a_crossSectionVector[i1] += factor1 * totalCrossSection1[i1] + factor2 * totalCrossSection2[i1];
    }
}

/* *********************************************************************************************************//**
 * Returns the requested reaction's multi-group cross section for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_reactionIndex       [in]    The index of the reaction.
 * @param a_hashIndex           [in]    The multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_sampling            [in]    If *true*, use augmented cross sections.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsMultiGroup::reactionCrossSection( int a_reactionIndex, int a_hashIndex, double a_temperature, bool a_sampling ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double cross_section;

    if( a_temperature <= m_temperatures[0] ) {
        cross_section = m_heatedCrossSections[0]->reactionCrossSection( a_reactionIndex, a_hashIndex, a_sampling ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        cross_section = m_heatedCrossSections.back( )->reactionCrossSection( a_reactionIndex, a_hashIndex, a_sampling ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        cross_section = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->reactionCrossSection( a_reactionIndex, a_hashIndex, a_sampling )
                            + fraction * m_heatedCrossSections[i1]->reactionCrossSection( a_reactionIndex, a_hashIndex, a_sampling );
    }

    return( cross_section );
}

/* *********************************************************************************************************//**
 * Returns the requested reaction's multi-group cross section for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_reactionIndex       [in]    The index of the reaction.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy_in           [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsMultiGroup::reactionCrossSection( int a_reactionIndex, double a_temperature, double a_energy_in ) const {

    int energyIndex = binarySearchVector( a_energy_in, m_projectileMultiGroupBoundariesCollapsed );

    if( energyIndex < 0 ) {
        energyIndex = 0;
        if( energyIndex == -1 ) energyIndex = static_cast<int>( m_projectileMultiGroupBoundariesCollapsed.size( ) ) - 2;
    }

    return( reactionCrossSection( a_reactionIndex, energyIndex, a_temperature, false ) );
}

/* *********************************************************************************************************//**
 * Returns the multi-group deposition energy for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_hashIndex           [in]    The multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 *
 * @return                              The deposition energy.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsMultiGroup::depositionEnergy( int a_hashIndex, double a_temperature ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double deposition_energy;

    if( a_temperature <= m_temperatures[0] ) {
        deposition_energy = m_heatedCrossSections[0]->depositionEnergy( a_hashIndex ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        deposition_energy = m_heatedCrossSections.back( )->depositionEnergy( a_hashIndex ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        deposition_energy = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->depositionEnergy( a_hashIndex )
                            + fraction * m_heatedCrossSections[i1]->depositionEnergy( a_hashIndex );
    }

    return( deposition_energy );
}

/* *********************************************************************************************************//**
 * Returns the multi-group deposition momentum for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_hashIndex           [in]    The multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 *
 * @return                              The deposition energy.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsMultiGroup::depositionMomentum( int a_hashIndex, double a_temperature ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double deposition_momentum;

    if( a_temperature <= m_temperatures[0] ) {
        deposition_momentum = m_heatedCrossSections[0]->depositionMomentum( a_hashIndex ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        deposition_momentum = m_heatedCrossSections.back( )->depositionMomentum( a_hashIndex ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );
        deposition_momentum = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->depositionMomentum( a_hashIndex )
                            + fraction * m_heatedCrossSections[i1]->depositionMomentum( a_hashIndex );
    }

    return( deposition_momentum );
}

/* *********************************************************************************************************//**
 * Returns the multi-group production energy for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_hashIndex           [in]    The multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 *
 * @return                              The deposition energy.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsMultiGroup::productionEnergy( int a_hashIndex, double a_temperature ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double production_energy;

    if( a_temperature <= m_temperatures[0] ) {
        production_energy = m_heatedCrossSections[0]->productionEnergy( a_hashIndex ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        production_energy = m_heatedCrossSections.back( )->productionEnergy( a_hashIndex ); }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );

        production_energy = ( 1. - fraction ) * m_heatedCrossSections[i1-1]->productionEnergy( a_hashIndex )
                            + fraction * m_heatedCrossSections[i1]->productionEnergy( a_hashIndex );
    }

    return( production_energy );
}

/* *********************************************************************************************************//**
 * Returns the multi-group gain for particle with index *a_particleIndex*. If no particle is found, a Vector of all 0's is returned.
 *
 * @param a_hashIndex           [in]    The multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_particleIndex       [in]    The index of the particle whose gain is to be returned.
 *
 * @return                              The multi-group gain.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsMultiGroup::gain( int a_hashIndex, double a_temperature, int a_particleIndex ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );

    if( a_temperature <= m_temperatures[0] ) {
        return( m_heatedCrossSections[0]->gain( a_particleIndex, a_hashIndex ) ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        return( m_heatedCrossSections.back( )->gain( a_particleIndex, a_hashIndex ) );
    }

    for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
    double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );

    double gain1 = m_heatedCrossSections[i1-1]->gain( a_particleIndex, a_hashIndex );
    double gain2 = m_heatedCrossSections[i1]->gain( a_particleIndex, a_hashIndex );

    return( ( 1. - fraction ) * gain1 + fraction * gain2 );
}

/* *********************************************************************************************************//**
 * Returns the multi-group gain for particle with intid *a_particleIntid*. If no particle is found, a Vector of all 0's is returned.
 *
 * @param a_hashIndex           [in]    The multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_particleIntid       [in]    The intid of the particle whose gain is to be returned.
 *
 * @return                              The multi-group gain.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double HeatedCrossSectionsMultiGroup::gainViaIntid( int a_hashIndex, double a_temperature, int a_particleIntid ) const {

    int i1, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );

    if( a_temperature <= m_temperatures[0] ) {
        return( m_heatedCrossSections[0]->gainViaIntid( a_particleIntid, a_hashIndex ) ); }
    else if( a_temperature >= m_temperatures.back( ) ) {
        return( m_heatedCrossSections.back( )->gainViaIntid( a_particleIntid, a_hashIndex ) );
    }

    for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
    double fraction = ( a_temperature - m_temperatures[i1-1] ) / ( m_temperatures[i1] - m_temperatures[i1-1] );

    double gain1 = m_heatedCrossSections[i1-1]->gainViaIntid( a_particleIntid, a_hashIndex );
    double gain2 = m_heatedCrossSections[i1]->gainViaIntid( a_particleIntid, a_hashIndex );

    return( ( 1. - fraction ) * gain1 + fraction * gain2 );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs index *a_particleIndex*.
 *
 * @param a_particleIndex       [in]    The PoPs index of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionsMultiGroup::setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {

    for( auto iter = m_heatedCrossSections.begin( ); iter != m_heatedCrossSections.end( ); ++iter ) (*iter)->setUserParticleIndex( a_particleIndex, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * Updates the m_userParticleIndex to *a_userParticleIndex* for all particles with PoPs intid *a_particleIntid*.
 *
 * @param a_particleIntid       [in]    The intid of the particle whose user index is to be set.
 * @param a_userParticleIndex   [in]    The particle index specified by the user.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionsMultiGroup::setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {

    for( auto iter = m_heatedCrossSections.begin( ); iter != m_heatedCrossSections.end( ); ++iter ) (*iter)->setUserParticleIndexViaIntid( a_particleIntid, a_userParticleIndex );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void HeatedCrossSectionsMultiGroup::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_VECTOR_DOUBLE( m_temperatures, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_thresholds, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_multiGroupThresholdIndex, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_projectileMultiGroupBoundariesCollapsed, a_buffer, a_mode );

    std::size_t vectorSize = m_heatedCrossSections.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_heatedCrossSections.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_heatedCrossSections.internalSize();
    for( std::size_t memberIndex = 0; memberIndex < vectorSize; ++memberIndex ) {
        if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
            if( a_buffer.m_placement != nullptr ) {
                m_heatedCrossSections[memberIndex] = new(a_buffer.m_placement) HeatedCrossSectionMultiGroup;
                a_buffer.incrementPlacement( sizeof( HeatedCrossSectionMultiGroup ) ); }
            else {
                m_heatedCrossSections[memberIndex] = new HeatedCrossSectionMultiGroup;
            }
        }
        if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof( HeatedCrossSectionMultiGroup ) );
        }
        m_heatedCrossSections[memberIndex]->serialize( a_buffer, a_mode );
    }
}

/* *********************************************************************************************************//**
 * This method writes the multi-group data at temperature index *a_temperatureIndex* to *a_file*.
 *
 * @param a_file                [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_temperatureIndex    [in]    The index of the temperature whose data are written.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionsMultiGroup::write( FILE *a_file, int a_temperatureIndex ) const {

    if( a_temperatureIndex < 0 ) return;

    std::size_t temperatureIndex = (std::size_t) a_temperatureIndex;
    if( temperatureIndex >= m_temperatures.size( ) ) return;

    printf( "HeatedCrossSectionsMultiGroup::write for temperature %.4e\n", m_temperatures[temperatureIndex] );

    fprintf( a_file, "    boundaries index                           " );
    std::string space( 14, ' ' );
    for( std::size_t index = 0; index < m_projectileMultiGroupBoundariesCollapsed.size( ); ++index ) {
        fprintf( a_file, "%s%6zu", space.c_str( ), index );
    }
    fprintf( a_file, "\n" );
    writeVector( a_file, "boundaries", 0, m_projectileMultiGroupBoundariesCollapsed );

    m_heatedCrossSections[temperatureIndex]->write( a_file );
}

/* *********************************************************************************************************//**
 * This method calls write for every temperature dataset in *this* with the file type stdout.
 ***********************************************************************************************************/

LUPI_HOST void HeatedCrossSectionsMultiGroup::print( ) const {

    for( std::size_t index = 0; index < m_heatedCrossSections.size( ); ++index ) write( stdout, index );
}

/* *********************************************************************************************************//**
 * Sets all elements of *a_vector* to 0.0 if *a_zeroReactions* is true, otherwise does nothing.
 *
 * @param a_vector              [in]    The vector to zero if *a_zeroReactions* is true.
 * @param a_zeroReactions       [in]    If true all elements of *a_vector* are set to 0.0.
 ***********************************************************************************************************/

static LUPI_HOST void checkZeroReaction( GIDI::Vector &a_vector, bool a_zeroReactions ) {

    if( a_zeroReactions ) a_vector.setToValueInFlatRange( 0, a_vector.size( ), 0.0 );
}

/* *********************************************************************************************************//**
 * Collapses the data in *a_vector* and calls *checkZeroReaction*.
 *
 * @param a_vector              [in]    The vector to collapse.
 * @param a_settings            [in]    Used to pass user options to the *this* to instruct it which data are desired.
 * @param a_particles           [in]    List of transporting particles and their information (e.g., multi-group boundaries and fluxes).
 * @param a_temperature         [in]    The temperature or the material.
 * @param a_zeroReactions       [in]    If true all elements of the returned **GIDI::Vector** are set to 0.0.
 ***********************************************************************************************************/

static LUPI_HOST GIDI::Vector collapseAndcheckZeroReaction( GIDI::Vector &a_vector, Transporting::MC const &a_settings, 
                GIDI::Transporting::Particles const &a_particles, LUPI_maybeUnused double a_temperature, bool a_zeroReactions ) {

    GIDI::Vector vector = GIDI::collapse( a_vector, a_settings, a_particles, 0.0 );
    checkZeroReaction( vector, a_zeroReactions );

    return( vector );
}

/* *********************************************************************************************************//**
 * @param a_file                [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_prefix              [in]    The prefix that starts the beginning of the line that the vector data are written to.
 * @param a_offset              [in]    Specifies the number of spaces to indent the line.
 * @param a_vector              [in]    The vector to write.
 ***********************************************************************************************************/

static void writeVector( FILE *a_file, std::string const &a_prefix, int a_offset, Vector<double> const &a_vector ) {

    std::string indent( 20 * a_offset, ' ' );

    std::string fmt = LUPI::Misc::argumentsToString( "    %%-%ds (%%4d) :: %%s", 40 );
    fprintf( a_file, fmt.c_str( ), a_prefix.c_str( ), (int) a_vector.size( ), indent.c_str( ) );
    for( Vector<double>::const_iterator iter = a_vector.begin( ); iter != a_vector.end( ); ++iter ) fprintf( a_file, " %19.11e", *iter );
    fprintf( a_file, "\n" );
}

}
