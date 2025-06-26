/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "MCGIDI.hpp"

#ifndef MCGIDI_CrossSectionLinearSubSearch
    #ifndef MCGIDI_CrossSectionBinarySubSearch
        #define MCGIDI_CrossSectionBinarySubSearch
    #endif
#endif

namespace MCGIDI {

namespace Sampling {

/* *********************************************************************************************************//**
 * This function returns the index in *a_energies* where *a_energy* lies between the returned index and the next index.
 * The returned index must lie between a_hashIndices[a_hashIndex] and a_hashIndices[a_hashIndex+1].
 * If *a_energy* is below the domain of *a_energies*, 0 is returned. If *a_energy* is above the domain of *a_energies*, 
 * the size of *a_energies* minus 2 is returned.
 * The argument *a_energyFraction* the weight for the energy at the returned index with the next index getting weighting 1 minus
 * *a_energyFraction*.
 *
 * @param a_hashIndex           [in]    The index in *a_hashIndices* where the index in *a_energy* must be bound by it and the next index in *a_hashIndex*.
 * @param a_hashIndices         [in]    The list of hash indices.
 * @param a_energies            [in]    The list of energies.
 * @param a_energy              [in]    The energy whose index is requested.
 * @param a_energyFraction      [in]    This represents the weighting to apply to the two bounding energies.
 *
 * @return                              The index bounding *a_energy* in the member *a_energies*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE int evaluationForHashIndex( int a_hashIndex, Vector<int> const &a_hashIndices, double a_energy, 
                Vector<double> const &a_energies, double *a_energyFraction ) {

    *a_energyFraction = 1.0;

    if( a_energy <= a_energies[0] ) return( 0 );
    if( a_energy >= a_energies.back( ) ) {
        *a_energyFraction = 0.0;
        return( (int) ( a_energies.size( ) - 2 ) );
    }

    int index1 = a_hashIndices[a_hashIndex];

#ifdef MCGIDI_CrossSectionLinearSubSearch
    while( a_energies[index1] > a_energy ) --index1;            // Make sure the calls gave the correct *a_hashIndex*.
    while( a_energies[index1] < a_energy ) ++index1;
    --index1;
#endif

#ifdef MCGIDI_CrossSectionBinarySubSearch
    int index2 = a_hashIndices[a_hashIndex];
    int index3 = (int) a_energies.size( ) - 1;
    if( ( a_hashIndex + 1 ) < (int) a_hashIndices.size( ) ) index3 = a_hashIndices[a_hashIndex+1] + 1;
    if( index3 == (int) a_energies.size( ) ) --index3;
    if( index2 != index3 ) index2 = binarySearchVectorBounded( a_energy, a_energies, index2, index3, false );
#endif

#ifdef MCGIDI_CrossSectionBinarySubSearch
    #ifdef MCGIDI_CrossSectionLinearSubSearch
        if( index1 != index2 ) {
            std::cerr << "Help " << index1 << "  " << index2 << std::endl;
        }
    #endif
    index1 = index2;
#endif

    *a_energyFraction = ( a_energies[index1+1] - a_energy ) / ( a_energies[index1+1] - a_energies[index1] );

    return( index1 );
}

namespace Upscatter {

/*! \class ModelDBRC_data
 * This class is used to store the cross section for the elastic scattering upscatter model B with Doppler Broadening
 * Rejection Correction (DBRC) with enum MCGIDI::Sampling::Upscatter::DBRC.
 */

/* *********************************************************************************************************//**
 *
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ModelDBRC_data::ModelDBRC_data( ) :
        m_neutronMass( 0.0 ),
        m_targetMass( ),
        m_energies( ),
        m_crossSections( ),
        m_hashIndices( ) {

}

/* *********************************************************************************************************//**
 *
 ***********************************************************************************************************/

LUPI_HOST ModelDBRC_data::ModelDBRC_data( double a_neutronMass, double a_targetMass, Vector<double> const  &a_energies, Vector<double> const &a_crossSections,
                DomainHash const &a_domainHash ) :
        m_neutronMass( a_neutronMass ),
        m_targetMass( a_targetMass ),
        m_energies( a_energies ),
        m_crossSections( a_crossSections ),
        m_hashIndices( a_domainHash.map( a_energies ) ),
        m_domainHash( 4000, 1e-8, 10 ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ModelDBRC_data::~ModelDBRC_data( ) {

}

/* *********************************************************************************************************//**
 * This method returns the cross section evaluate at the projectile speed *a_speed*.
 *
 * @param a_temperature         [in]    The temperature of the target.
 *
 * @return                              The thermal speed of the target.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ModelDBRC_data::evaluate( double a_energy ) {

    double energyFraction;
    int hashIndex = m_domainHash.index( a_energy );
    int index = evaluationForHashIndex( hashIndex, m_hashIndices, a_energy, m_energies, &energyFraction );

    return( energyFraction * m_crossSections[index] + ( 1.0 - energyFraction ) * m_crossSections[index+1] );
}

/* *********************************************************************************************************//**
 * This method returns the thermal speed of the target with temperature *a_temperature*.
 *
 * @param a_temperature         [in]    The temperature of the target.
 *
 * @return                              The thermal speed of the target.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ModelDBRC_data::targetThermalSpeed( double a_temperature ) {

    return( sqrt( 2.0 * a_temperature / m_targetMass ) );
}

/* *********************************************************************************************************//**
 * This method sets *a_crossSectionMin* and *a_crossSectionMax* to the minimum and maximum cross section values
 * for the cross section in the window (*a_speed* - 4 * *a_targetThermalSpeed*) < speed < (*a_speed* + 4 * *a_targetThermalSpeed*).
 * If (a_speed - 4 * a_targetThermalSpeed) is less than the lowest speed in the data, then it is replaced with the
 * lowest speed.
 *
 * @param a_energy              [in]    The energy of the projectile (i.e., incident neutron).
 * @param a_targetThermalSpeed  [in]    The thermal speed of the target.
 *
 * @return                              The maximum cross section in the search window.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ModelDBRC_data::crossSectionMax( double a_energy, double a_targetThermalSpeed ) {

    double crossSectionMax2 = 0.0;
    double energyFraction;
    double a_speed = MCGIDI_particleBeta( m_neutronMass, a_energy );

    double speedMin = a_speed - 4 * a_targetThermalSpeed;
    if( speedMin < 0.0 ) speedMin = 0.0;
    double energyMin = 0.5 * m_neutronMass * speedMin * speedMin;
    int hashIndex = m_domainHash.index( energyMin );
    int indexMin = evaluationForHashIndex( hashIndex, m_hashIndices, energyMin, m_energies, &energyFraction );

    double speedMax = a_speed + 4 * a_targetThermalSpeed;
    double energyMax = 0.5 * m_neutronMass * speedMax * speedMax;
    hashIndex = m_domainHash.index( energyMax );
    int indexMax = evaluationForHashIndex( hashIndex, m_hashIndices, energyMax, m_energies, &energyFraction );
    if( indexMax < static_cast<int>( m_energies.size( ) ) ) ++indexMax;

    for( int index = indexMin; index < indexMax; ++index ) {
        if( crossSectionMax2 < m_crossSections[index] ) crossSectionMax2 = m_crossSections[index];
    }

    return( crossSectionMax2 );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ModelDBRC_data::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_DOUBLE( m_neutronMass, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_targetMass, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_energies, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_crossSections, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_INT( m_hashIndices, a_buffer, a_mode );
    m_domainHash.serialize( a_buffer, a_mode );
}

/* *********************************************************************************************************//**
 * This method serializes data for a *ModelDBRC_data* instance.
 *
 * @param a_modelDBRC_data      [in/out]    A pointer to the **ModelDBRC_data** instance to serialize.
 * @param a_buffer              [in]        The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]        Specifies the action of this method.
 *
 * @returns                                 A pointer the serialized **ModelDBRC_data** instance.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ModelDBRC_data *serializeModelDBRC_data( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, ModelDBRC_data *a_modelDBRC_data ) {

    bool haveDBRC = a_modelDBRC_data != nullptr;
    DATA_MEMBER_CAST( haveDBRC, a_buffer, a_mode, bool );

    if( haveDBRC ) {
        if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
            if (a_buffer.m_placement != nullptr) {
                a_modelDBRC_data = new(a_buffer.m_placement) ModelDBRC_data;
                a_buffer.incrementPlacement( sizeof( ModelDBRC_data ) ); }
            else {
                a_modelDBRC_data = new ModelDBRC_data;
            }
        }

        if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
            a_buffer.incrementPlacement( sizeof( ModelDBRC_data ) );
        }

        a_modelDBRC_data->serialize( a_buffer, a_mode );
    }

    return( a_modelDBRC_data );
}

}

/*
=========================================================
*/
LUPI_HOST_DEVICE ClientRandomNumberGenerator::ClientRandomNumberGenerator( double (*a_generator)( void * ), void *a_state ) :
        m_generator( a_generator ),
        m_state( a_state ) {
}

/*
=========================================================
*/
LUPI_HOST_DEVICE ClientCodeRNGData::ClientCodeRNGData( double (*a_generator)( void * ), void *a_state ) :
        ClientRandomNumberGenerator( a_generator, a_state ) {
}

/*
=========================================================
*/
LUPI_HOST_DEVICE Input::Input( bool a_wantVelocity, Upscatter::Model a_upscatterModel ) :
        m_wantVelocity( a_wantVelocity ),
        m_upscatterModel( a_upscatterModel ) {

}

}           // End of namespace Sampling.

}           // End of namespace MCGIDI.
