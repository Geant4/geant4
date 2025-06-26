/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <climits>
#include <PoPI.hpp>

namespace PoPI {

/* *********************************************************************************************************//**
 * Returns an integer representing the particle's famuly *a_family*.
 *
 * @param a_isAnti          [in]    If **true** particle is an anti-particle and otherwise its a particle.
 ***********************************************************************************************************/

int family2Integer( Particle_class a_family ) {

    if( a_family == Particle_class::nucleus ) return( -1 );
    if( a_family == Particle_class::nuclide ) return( -2 );
    if( a_family == Particle_class::gaugeBoson ) return( 0 );
    if( a_family == Particle_class::lepton ) return( 1 );
    if( a_family == Particle_class::baryon ) return( 2 );
    if( a_family == Particle_class::nuclideMetaStable ) return( 50 );
    if( a_family == Particle_class::nucleusMetaStable ) return( 60 );
    if( a_family == Particle_class::ENDL_fissionProduct ) return( 99 );

    return( -3 );
}

/* *********************************************************************************************************//**
 * This function is for internal use.
 * Returns the intid for the particle of family *a_family* with family indentifier *a_SSSSSSS*. If the return value is -1
 * the family is not supported by this function.
 *
 * @param a_isAnti          [in]    If **true** particle is an anti-particle and otherwise its a particle.
 * @param a_family          [in]    The particle's family.
 * @param a_SSSSSSS         [in]    The particle's indentifier within its family.
 *
 * @return                          The intid for the particle.
 ***********************************************************************************************************/

int intidHelper( bool a_isAnti, Particle_class a_family, int a_SSSSSSS ) {

    int sign = a_isAnti ? -1 : 1;

    int intid = family2Integer( a_family );
    if( intid < 0 ) return( -1 );
    intid += 100;
    intid *= 10000000;

    return( sign * ( intid + a_SSSSSSS ) );
}

/*! \class ParseIntidInfo
 * This class represents **PoPs** nucleus instance.
 */

/* *********************************************************************************************************//**
 * Constructor that parses *a_intid* into its components and sets members per *a_intid*. If *m_III*, *m_ZZZ*, *m_AAA* and
 * *m_metaStableIndex* are positive (greater than or equal to 0), then the particle is a nuclear particle and 
 * *m_isNuclear* is *true*, otherwise the particle is not a nuclear particle and *m_isNuclear* is *false*.
 * Note, even if *m_metaStableIndex* > 0 (i.e., particle is a nuclear meta-stable), *m_III* is as expected. For
 * example, for intid = 481095242, *m_metaStableIndex* is 1 and *m_III* is 481.
 *
 * If *m_family* is **Particle_class::unknown** then all other members are undefined.
 *
 * @param a_intid           [in]    The intid for the particle to parse.
 * @param a_GRIN_mode       [in]    For the GRIN project, nuclear levels go beyond 499, so this flag causes III >= 500 to be treated as a nuclide.
 ***********************************************************************************************************/

ParseIntidInfo::ParseIntidInfo( int a_intid, bool a_GRIN_mode ) :
        m_intid( a_intid ),
        m_family( Particle_class::unknown ),
        m_isAnti( a_intid < 0 ),
        m_isNuclear( false ),
        m_AAA( -1 ),
        m_ZZZ( -1 ),
        m_III( -1 ),
        m_metaStableIndex( -1 ),
        m_generation( -1 ),
        m_isNeutrino( false ),
        m_baryonGroup( -1 ),
        m_baryonId( -1 ),
        m_familyId( -1 ) {

    int intidAbs = std::abs( a_intid );

    bool nuclearLike = intidAbs / 1000000000 == 0;
    int family = (intidAbs / 10000000) % 100;
    int SSSSSSS = intidAbs % 10000000;

    if(      nuclearLike ) {
        m_AAA = intidAbs % 1000;
        m_ZZZ = intidAbs % 1000000 / 1000;
        m_III = intidAbs % 1000000000 / 1000000;

        m_nuclearLevelIndex = m_III;
        if( ( m_III < 500 ) || a_GRIN_mode ) {
            m_family = Particle_class::nuclide; }
        else {
            m_nuclearLevelIndex -= 500;
            m_family =  Particle_class::nucleus;
        } }
    else {
        int topFamilyDigid = family / 10;
        if( ( topFamilyDigid == 5 ) || ( topFamilyDigid == 6 ) ) {
            m_family = topFamilyDigid == 5 ? Particle_class::nuclideMetaStable : Particle_class::nucleusMetaStable;
            m_AAA = intidAbs % 1000;
            m_ZZZ = intidAbs % 1000000 / 1000;
            m_metaStableIndex = intidAbs % 100000000 / 1000000; }
        else {
            m_familyId = SSSSSSS;
            if(      family == 0 ) {
                m_family =  Particle_class::gaugeBoson; }
            else if( family == 1 ) {
                int neutronoFlag = ( SSSSSSS % 100 ) / 10;
                if( neutronoFlag > 1 ) return;                      // Invalid particle.

                m_family =  Particle_class::lepton;
                m_generation = SSSSSSS % 10;
                m_isNeutrino = neutronoFlag != 0; }
            else if( family == 2 ) {
                m_family =  Particle_class::baryon;
                m_baryonGroup = SSSSSSS / 1000000;
                m_baryonId = SSSSSSS % 1000000; }
            else if( family == 98 ) {
                m_family =  Particle_class::TNSL; }
            else if( family == 99 ) {
                m_family =  Particle_class::ENDL_fissionProduct;
            }
        }
    }
}

/* *********************************************************************************************************//**
 * Returns the GNDS PoPs id for *this*. If particles is unknown, an empty string is returned.
 *
 ***********************************************************************************************************/

std::string ParseIntidInfo::id( ) {

    std::string pid;

    if( ( m_family == Particle_class::nuclide ) || ( m_family == Particle_class::nucleus ) || ( m_family == Particle_class::nuclideMetaStable ) ||
                ( m_family == Particle_class::nucleusMetaStable ) ) {
        bool isNucleus = ( m_family == Particle_class::nucleus ) || ( m_family == Particle_class::nucleusMetaStable );
        pid = chemicalElementInfoFromZ( m_ZZZ, true, isNucleus );
        if( pid != "" ) {
            pid += LUPI::Misc::argumentsToString( "%d", m_AAA );

            if( ( m_family == Particle_class::nuclideMetaStable ) || ( m_family == Particle_class::nucleusMetaStable ) ) {
                    pid += LUPI::Misc::argumentsToString( "_m%d", m_metaStableIndex ); }
            else {
                int III = m_III;
                if( m_family == Particle_class::nucleus ) III -= 500;
                if( III != 0 ) pid += LUPI::Misc::argumentsToString( "_e%d", III );
            }
        } }
    else if( m_family == Particle_class::gaugeBoson ) {
        if( m_familyId == 0 ) pid = IDs::photon; }
    else if( m_family == Particle_class::lepton ) {
        if( m_generation == 0 ) {
            if( !m_isNeutrino ) pid = IDs::electron;
        } }
    else if( m_family == Particle_class::baryon ) {
        if( m_baryonGroup == 0 ) {
            switch( m_baryonId ) {
            case 0:
                pid = IDs::neutron;
                break;
            case 1:
                pid = IDs::proton;
                break;
            default:
                break;
            }
        } }
    else if( m_family == Particle_class::ENDL_fissionProduct ) {
        if( m_familyId == 99120 ) {
            pid = IDs::FissionProductENDL99120; }
        else if( m_familyId == 99125 ) {
            pid = IDs::FissionProductENDL99125;
        }
    }

    if( ( pid.size( ) > 0 ) &&  m_isAnti ) pid += IDs::anti;

    return( pid );
}

}
