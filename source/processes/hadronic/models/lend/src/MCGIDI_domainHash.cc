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

/*! \class DomainHash
 * This class stores the data needed for logarithmic hash look up of a domain. This is used to find a cross section given a projectile's energy.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE DomainHash::DomainHash( ) :
        m_bins( 0 ),
        m_domainMin( 0.0 ),
        m_domainMax( 0.0 ),
        m_u_domainMin( 0.0 ),
        m_u_domainMax( 0.0 ),
        m_inverse_du( 0.0 ) {

}

/* *********************************************************************************************************//**
 * @param a_bins                [in]    The number of bins for the hahs function.
 * @param a_domainMin           [in]    The minimum value of the energy domain for the hash function.
 * @param a_domainMax           [in]    The maximum value of the energy domain for the hash function.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE DomainHash::DomainHash( std::size_t a_bins, double a_domainMin, double a_domainMax ) :
        m_bins( a_bins ),
        m_domainMin( a_domainMin ),
        m_domainMax( a_domainMax ),
        m_u_domainMin( log( a_domainMin ) ),
        m_u_domainMax( log( a_domainMax ) ),
        m_inverse_du( static_cast<double>( a_bins ) / ( m_u_domainMax - m_u_domainMin ) ) {

}

/* *********************************************************************************************************//**
 * @param a_domainHash          [in]    The DomainHash instance to copy.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE DomainHash::DomainHash( DomainHash const &a_domainHash ) :
        m_bins( a_domainHash.bins( ) ),
        m_domainMin( a_domainHash.domainMin( ) ),
        m_domainMax( a_domainHash.domainMax( ) ), 
        m_u_domainMin( a_domainHash.u_domainMin( ) ),
        m_u_domainMax( a_domainHash.u_domainMax( ) ),
        m_inverse_du( a_domainHash.inverse_du( ) ) {
}

/* *********************************************************************************************************//**
 * This method returns the hash index given the domain value *a_domain*. If *a_domain* is less than *m_domainMin*,
 * the returned index is 0. If *a_domain* is greater than *m_domainMax*, the returned index is *m_bins* + 1.
 * Otherwise, the returned index is in the range [1, *m_bins*].
 *
 * @param a_domain              [in]    The domain value that the hash index is to be returned for.
 *
 * @return                              The hash index.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE std::size_t DomainHash::index( double a_domain ) const {

    if( a_domain < m_domainMin ) return( 0 );
    if( a_domain > m_domainMax ) return( m_bins + 1 );
    double dIndex = m_inverse_du * ( log( a_domain ) - m_u_domainMin ) + 1;
    return( static_cast<std::size_t>( dIndex ) );
}

/* *********************************************************************************************************//**
 * This method returns the hash indices for the requested domain values *a_domainValues*.
 *
 * @param a_domainValues        [in]    The domain values.
 *
 * @return                              The hash indices.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Vector<std::size_t> DomainHash::map( Vector<double> const &a_domainValues ) const {

    std::size_t i1, size( a_domainValues.size( ) );
    Vector<std::size_t> indices( m_bins + 2, static_cast<std::size_t>( 0 ) );
    std::size_t lastIndex = 0, currentIndex, i2 = 1;

    for( i1 = 0; i1 < size; ++i1 ) {
        currentIndex = index( a_domainValues[i1] );
        if( currentIndex != lastIndex ) {
            for( ; lastIndex < currentIndex; ++lastIndex, ++i2 ) {
                indices[i2] = i1 - 1;
                if( i1 == 0 ) indices[i2] = 0;          // Special case.
            }
        }
    }
    for( ; i2 < ( m_bins + 2 ); ++i2 ) indices[i2] = indices[i2-1];
    return( indices );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void DomainHash::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_SIZE_T( m_bins, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_domainMin, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_domainMax, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_u_domainMin, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_u_domainMax, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_inverse_du, a_buffer, a_mode );
}

/* *********************************************************************************************************//**
 * Prints the contents of *this*.
 *
 * @param a_printValues         [in]    If true, the domain values that divide the hash indices are also printed.
 ***********************************************************************************************************/

LUPI_HOST void DomainHash::print( bool a_printValues ) const {
#ifndef __CUDA_ARCH__
    std::cout << "bins = " << m_bins << std::endl;
    std::cout << "    m_domainMin = " << m_domainMin << "  << m_domainMax = " << m_domainMax << std::endl;
    std::cout << "    m_u_domainMin = " << m_u_domainMin << "  << m_u_domainMax = " << m_u_domainMax << std::endl;
    std::cout << "    m_inverse_du = " << m_inverse_du << std::endl;
    if( a_printValues ) {
        double domain = m_domainMin, factor = pow( m_domainMax / m_domainMin, 1. / static_cast<double>( m_bins ) );

        for( std::size_t i1 = 0; i1 < bins( ); ++i1, domain *= factor ) {
            std::cout << LUPI::Misc::argumentsToString( " %14.7e", domain );
            if( ( ( i1 + 1 ) % 10 ) == 0 ) std::cout << std::endl;
        }
        std::cout << LUPI::Misc::argumentsToString( " %14.7e", m_domainMax ) << std::endl;
    }
#endif
}

/*! \class MultiGroupHash
 * This class stores a multi-group boundaries and has a method *index* that returns an index of the group for the requested domain value.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE MultiGroupHash::MultiGroupHash( ) {

}

/* *********************************************************************************************************//**
 * @param a_boundaries          [in]    The list of multi-group boundaries.
 ***********************************************************************************************************/

LUPI_HOST MultiGroupHash::MultiGroupHash( std::vector<double> a_boundaries ) :
        m_boundaries( a_boundaries ) {

}

/* *********************************************************************************************************//**
 * @param a_boundaries          [in]    The list of multi-group boundaries.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE MultiGroupHash::MultiGroupHash( Vector<double> a_boundaries ) :
        m_boundaries( a_boundaries ) {

}

/* *********************************************************************************************************//**
 * This constructor gets the list of multi-group boundaries from the first GIDI::Styles::MultiGroup of *a_protare*.
 * It calls MultiGroupHash::initialize to set up *this*.
 *
 * @param a_protare             [in]    The GIDI::Protare containing the GIDI::Styles::MultiGroup style.
 * @param a_temperatureInfo     [in]    This is used to determine the multi-group boundaries.
 * @param a_particleID          [in]    The PoPs' id of the particle whose multi-group boundaries are desired.
 ***********************************************************************************************************/

LUPI_HOST MultiGroupHash::MultiGroupHash( GIDI::Protare const &a_protare, GIDI::Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_particleID ) {

    initialize( a_protare, a_temperatureInfo, a_particleID );
}

/* *********************************************************************************************************//**
 * This constructor gets the list of multi-group boundaries from the GIDI::Particle of *a_particles* that is the projectile.
 *
 * @param a_protare             [in]    The GIDI::Protare containing the GIDI::Styles::MultiGroup style.
 * @param a_particles           [in]    The list of transportable particles.
 ***********************************************************************************************************/

LUPI_HOST MultiGroupHash::MultiGroupHash( GIDI::Protare const &a_protare, GIDI::Transporting::Particles const &a_particles ) {

    GIDI::Transporting::Particle const &particle = *a_particles.particle( a_protare.projectile( ).pid( ) );

    m_boundaries = particle.multiGroup( ).boundaries( );
}

/* *********************************************************************************************************//**
 * This constructor gets the list of multi-group boundaries from the GIDI::Particle of *a_particles* that is the projectile.
 *
 * @param a_protare             [in]    The GIDI::Protare containing the GIDI::Styles::MultiGroup style.
 * @param a_particles           [in]    The list of transportable particles.
 ***********************************************************************************************************/

LUPI_HOST MultiGroupHash::MultiGroupHash( MultiGroupHash const &a_multiGroupHash ) :
        m_boundaries( a_multiGroupHash.boundaries( ) ) {

}

/* *********************************************************************************************************//**
 * This method is used by several constructors to get the multi-group data.
 *
 * @param a_protare             [in]    The GIDI::Protare containing the GIDI::Styles::MultiGroup style.
 * @param a_temperatureInfo     [in]    This is used to determine the multi-group boundaries.
 * @param a_particleID          [in]    The PoPs' id of the particle whose multi-group boundaries are desired.
 ***********************************************************************************************************/

LUPI_HOST void MultiGroupHash::initialize( GIDI::Protare const &a_protare, GIDI::Styles::TemperatureInfo const &a_temperatureInfo, std::string a_particleID ) {

    if( a_particleID == "" ) a_particleID = a_protare.projectile( ).ID( );

    GIDI::Styles::Suite const &stylesSuite( a_protare.styles( ) );
    GIDI::Styles::HeatedMultiGroup const *heatedMultiGroupStyle1 = stylesSuite.get<GIDI::Styles::HeatedMultiGroup>( a_temperatureInfo.heatedMultiGroup( ) );

    m_boundaries = heatedMultiGroupStyle1->groupBoundaries( a_particleID );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void MultiGroupHash::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_VECTOR_DOUBLE( m_boundaries, a_buffer, a_mode );
}

}
