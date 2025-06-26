/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <limits.h>
#include <sstream>
#include <stdexcept>

#include "PoPI.hpp"

namespace PoPI {

std::string const IDs::photon = "photon";
std::string const IDs::electron = "e-";
std::string const IDs::neutron = "n";
std::string const IDs::proton = "p";
std::string const IDs::familiarPhoton = "g";
std::string const IDs::familiarDeuteron = "d";
std::string const IDs::familiarTriton = "t";
std::string const IDs::familiarHelion = "h";
std::string const IDs::familiarAlpha = "a";
std::string const IDs::FissionProductENDL99120 = "FissionProductENDL99120";
std::string const IDs::FissionProductENDL99125 = "FissionProductENDL99125";
std::string const IDs::anti = "_anti";

/* *********************************************************************************************************//**
 * Returns true if a_formatVersion is a format supported by **PoPI** and false otherwise;
 *
 * @param a_formatVersion               [in]    The format version to check if it is supported.
 *  
 * @return                                      **true** if format is supported by **PoPI** and **false** otherwise.
 ***********************************************************************************************************/

bool supportedFormat( LUPI::FormatVersion const &a_formatVersion ) {

    if( a_formatVersion.format( ) == PoPI_formatVersion_0_1_Chars ) return( true );

    return( a_formatVersion.supported( ) );
}

/* *********************************************************************************************************//**
 * Added the XML end tag (e.g., "</tag>") with tag name *a_label* to the last srd::string in *a_XMLList*.
 *
 * @param a_XMLList                     [in]    The list whose last item is emended with an XML end tag.
 * @param a_label                       [in]    The name of the end tag.
 ***********************************************************************************************************/

void appendXMLEnd( std::vector<std::string> &a_XMLList, std::string const &a_label ) {

    std::string theEnd = "</" + a_label + ">";
    std::vector<std::string>::iterator iter = a_XMLList.end( );
    --iter;
    *iter += theEnd;
}

/* *********************************************************************************************************//**
 * This function returns the "special particle id" for the specified particle id (*a_id*). The id returned depends
 * on the *a_mode* argument. Currently, for all *a_id*'s but for those listed in the following table, * *a_id* is 
 * returned. For an *a_id* in the following table, find the column that contains the specfied *a_id*, the id returned
 * will be the id in that column whose row matches the specified *a_mode*.
 * 
   | a_mode       | ids -> |||||
   | ------------ | ----- | ----- | ----- | ----- | ---- |
   | familiar     | p     | d     | t     | h     | a    |
   | nuclide      | H1    | H2    | H3    | He3   | He4  |
   | nucleus      | h1    | h2    | h3    | he3   | he4  |
 *
 * @param a_mode                [in]    The mode which determines the returned id.
 * @param a_id                  [in]    The specified particle id.
 *
 * @return                              The special particle id.
 ***********************************************************************************************************/

std::string specialParticleID( SpecialParticleID_mode a_mode, std::string const &a_id ) {

    static std::string firstChars( "pdthaH" );
    int iid = 0;

    if( a_id.size( ) > 3 ) return( a_id );

    std::size_t index = firstChars.find( a_id[0] );

    if( index == std::string::npos ) return( a_id );

    if( a_id == "H" ) return( a_id );                    // Special case to avoid test of id1[0] two lines later.
    std::string id1( a_id );
    if( id1[0] == 'H' ) id1[0] = 'h';

    if(      id1 == IDs::proton || id1 == "h1" ) {
        iid = 1; }
    else if( id1 == IDs::familiarDeuteron || id1 == "h2" ) {
        iid = 2; }
    else if( id1 == IDs::familiarTriton   || id1 == "h3" ) {
        iid = 3; }
    else if( id1 == IDs::familiarHelion   || id1 == "he3" ) {
        iid = 4; }
    else if( id1 == IDs::familiarAlpha    || id1 == "he4" ) {
        iid = 5;
    }
    if( iid == 0 ) return( a_id );

    if( a_mode == SpecialParticleID_mode::familiar ) {
        return( firstChars.substr( iid-1, 1 ) ); }
    else {
        if(      iid == 1 ) {
            id1 = "h1"; }
        else if( iid == 2 ) {
            id1 = "h2"; }
        else if( iid == 3 ) {
            id1 = "h3"; }
        else if( iid == 4 ) {
            id1 = "he3"; }
        else {
            id1 = "he4";
        }
        if( a_mode == SpecialParticleID_mode::nuclide ) id1[0] = 'H';
    }
    return( id1 );
}

/* *********************************************************************************************************//**
 * Compares two particle ids and returns if they are the same particle. This methods using the name returned by the function
 * **specialParticleID** with the same SpecialParticleID_mode for each particle id. Ergo, "H1" is the same as "H1", "p" or "h1".
 *
 * @return                                  **true** if particles are the same and **false** otherwise.
 ***********************************************************************************************************/

bool compareSpecialParticleIDs( std::string const &a_id1, std::string const &a_id2 ) {

    return( specialParticleID( SpecialParticleID_mode::familiar, a_id1 ) == specialParticleID( SpecialParticleID_mode::familiar, a_id2 ) );
}

/* *********************************************************************************************************//**
 * Returns the Z (i.e., the atomic number) for a nuclear type particle; otherwise, 0 is returned. Currently, non-0
 * values are returned if *a_particle* is an isotope, nuclide or nucleus, or if it is a proton and *a_isNeutronProtonANucleon*
 * is **true**. Note, this is not the charge of the particle but its atomic number. For example, the atomic number
 * for an electron is 0 as it is not a nuclear type particle.
 *  
 * @param a_particle                    [in]    The PoPI::Base instance whose Z is returned.
 * @param a_isNeutronProtonANucleon     [in]    If **true** a proton is treated as a nucleus.
 *
 * @return                                  The Z of the particle
 ***********************************************************************************************************/

int particleZ( Base const &a_particle, bool a_isNeutronProtonANucleon ) {

    int Z = 0;

    if( a_particle.ID( ) == IDs::proton ) {
        if( a_isNeutronProtonANucleon ) Z = 1; }
    else if( a_particle.isNuclide( ) ) {
        Nuclide const &particle = (Nuclide const &) a_particle;
        Z = particle.Z( ); }
    else if( a_particle.isNucleus( ) ) {
        Nucleus const &particle = (Nucleus const &) a_particle;
        Z = particle.Z( ); }
    else if( a_particle.isChemicalElement( ) ) {
        ChemicalElement const &object = (ChemicalElement const &) a_particle;
        Z = object.Z( ); }
    else if( a_particle.isIsotope( ) ) {
        Isotope const &object = (Isotope const &) a_particle;
        Z = object.Z( );
    }

    return( Z );
}

/* *********************************************************************************************************//**
 * Uses the index *a_index* to look up the particle in *a_pops* and calls **particleZ** for that particle.
 *
 * @param a_pops                        [in]    The PoPs database to look up the particle.
 * @param a_index                       [in]    The index of the particle in *a_pops* whose Z value is returned.
 * @param a_isNeutronProtonANucleon     [in]    If **true** a proton is treated as a nucleus.
 *
 * @return                                      The Z returned by **particleZ( Base const &, bool )**.
 ***********************************************************************************************************/

int particleZ( Database const &a_pops, int a_index, bool a_isNeutronProtonANucleon ) {

    int Z = 0;
    Base const &base( a_pops.get<Base>( a_pops.final( a_index ) ) );

    if( base.isChemicalElement( ) ) {
        SymbolBase const &object2( static_cast<SymbolBase const &>( base ) );
        Z = particleZ( object2 ); }
    else {
        Particle const &particle( static_cast<Particle const &>( base ) );
        Z = particleZ( particle, a_isNeutronProtonANucleon );
    }

    return( Z );
}

/* *********************************************************************************************************//**
 * Uses the id *a_id* to look up the particle in *a_pops* and calls **particleZ** for that particle.
 *
 * @param a_pops                        [in]    The PoPs database to look up the particle.
 * @param a_id                          [in]    The id of the particle in *a_pops* whose Z value is returned.
 * @param a_isNeutronProtonANucleon     [in]    If **true** a proton is treated as a nucleus.
 *
 * @return                                      The Z returned by **particleZ( Base const &, bool )**.
 ***********************************************************************************************************/

int particleZ( Database const &a_pops, std::string const &a_id, bool a_isNeutronProtonANucleon ) {

    Base const &object( a_pops.get<Particle>( a_pops.final( a_id ) ) );

    return( particleZ( object, a_isNeutronProtonANucleon ) );
}

/* *********************************************************************************************************//**
 * Returns the A (i.e., atomic mass number) for a nuclear type particle; otherwise, 0 is returned. Currently, non-0
 * values are returned if *a_particle* is an isotope, nuclide or nucleus, or if it is a neutron or a proton and 
 * *a_isNeutronProtonANucleon* is **true**.
 *  
 * @param a_particle                    [in]    The PoPI::Base instance whose Z is returned.
 * @param a_isNeutronProtonANucleon     [in]    If **true** a proton is treated as a nucleus.
 *
 * @return                                  The Z of the particle
 ***********************************************************************************************************/

int particleA( Base const &a_particle, bool a_isNeutronProtonANucleon ) {

    int A = 0;

    if( a_particle.ID( ) == IDs::neutron ) {
        if( a_isNeutronProtonANucleon ) A = 1; }
    else if( a_particle.ID( ) == IDs::proton ) {
        if( a_isNeutronProtonANucleon ) A = 1; }
    else if( a_particle.isNuclide( ) ) {
        Nuclide const &particle = (Nuclide const &) a_particle;
        A = particle.A( ); }
    else if( a_particle.isNucleus( ) ) {
        Nucleus const &particle = (Nucleus const &) a_particle;
        A = particle.A( ); }
    else if( a_particle.isIsotope( ) ) {
        Isotope const &object = (Isotope const &) a_particle;
        A = object.A( );
    }

    return( A );
}

/* *********************************************************************************************************//**
 * Uses the index *a_index* to look up the particle in *a_pops* and calls **particleA** for that particle.
 *
 * @param a_pops                        [in]    The PoPs database to look up the particle.
 * @param a_index                       [in]    The index of the particle in *a_pops* whose A value is returned.
 * @param a_isNeutronProtonANucleon     [in]    If **true** a proton is treated as a nucleus.
 *
 * @return                                      The A returned by **particleA( Base const &, bool )**.
 ***********************************************************************************************************/

int particleA( Database const &a_pops, int a_index, bool a_isNeutronProtonANucleon ) {

    Base const &particle( a_pops.get<Base>( a_pops.final( a_index ) ) );

    return( particleA( particle, a_isNeutronProtonANucleon ) );
}

/* *********************************************************************************************************//**
 * Uses the id *a_id* to look up the particle in *a_pops* and calls **particleA** for that particle.
 *
 * @param a_pops                        [in]    The PoPs database to look up the particle.
 * @param a_id                          [in]    The id of the particle in *a_pops* whose A value is returned.
 * @param a_isNeutronProtonANucleon     [in]    If **true** a proton is treated as a nucleus.
 *
 * @return                                      The A returned by **particleA( Base const &, bool )**.
 ***********************************************************************************************************/

int particleA( Database const &a_pops, std::string const &a_id, bool a_isNeutronProtonANucleon ) {

    Base const &particle( a_pops.get<Base>( a_pops.final( a_id ) ) );

    return( particleA( particle, a_isNeutronProtonANucleon ) );
}

/* *********************************************************************************************************//**
 * Returns the ZA (i.e., 1000 * Z + A) for a nuclear type particle; otherwise, 0 is returned. Currently, non-0
 * values are returned if *a_particle* is an isotope, nuclide or nucleus, or if it is a neutron or a proton 
 * and *a_isNeutronProtonANucleon* is **true**.
 *  
 * @param a_particle                    [in]    The PoPI::Base instance whose Z is returned.
 * @param a_isNeutronProtonANucleon     [in]    If **true** a proton is treated as a nucleus.
 *
 * @return                                  The Z of the particle
 ***********************************************************************************************************/

int particleZA( Base const &a_particle, bool a_isNeutronProtonANucleon ) {

    int ZA = 0;

    if( a_particle.ID( ) == IDs::neutron ) {
        if( a_isNeutronProtonANucleon ) ZA = 1; }
    else {
        if( !a_particle.isChemicalElement( ) ) ZA = 1000 * particleZ( a_particle, a_isNeutronProtonANucleon ) + particleA( a_particle, a_isNeutronProtonANucleon );
    }

    return( ZA );
}

/* *********************************************************************************************************//**
 * Uses the index *a_index* to look up the particle in *a_pops* and calls **particleZA** for that particle.
 *
 * @param a_pops                        [in]    The PoPs database to look up the particle.
 * @param a_index                       [in]    The index of the particle in *a_pops* whose ZA value is returned.
 * @param a_isNeutronProtonANucleon     [in]    If **true** a proton is treated as a nucleus.
 *
 * @return                                      The ZA returned by **particleZA( Base const &, bool )**.
 ***********************************************************************************************************/

int particleZA( Database const &a_pops, int a_index, bool a_isNeutronProtonANucleon ) {

    Base const &particle( a_pops.get<Base>( a_pops.final( a_index ) ) );
    
    return( particleZA( particle, a_isNeutronProtonANucleon ) );
}

/* *********************************************************************************************************//**
 * Uses the id *a_id* to look up the particle in *a_pops* and calls **particleZA** for that particle.
 *
 * @param a_pops                        [in]    The PoPs database to look up the particle.
 * @param a_id                          [in]    The id of the particle in *a_pops* whose ZA value is returned.
 * @param a_isNeutronProtonANucleon     [in]    If **true** a proton is treated as a nucleus.
 *
 * @return                                      The ZA returned by **particleZA( Base const &, bool )**.
 ***********************************************************************************************************/

int particleZA( Database const &a_pops, std::string const &a_id, bool a_isNeutronProtonANucleon ) {

    Base const &particle( a_pops.get<Base>( a_pops.final( a_id ) ) );

    return( particleZA( particle, a_isNeutronProtonANucleon ) );
}

/* *********************************************************************************************************//**
 * Returns the meta-stable index if *a_particle* is a PoPI::MetaStable alias; otherwise, 0 is returned. 
 *
 * @param a_particle                    [in]    The PoPI::Base instance whose meta-stable index is returned.
 *
 * @return                                      The meta-stable index of the particle
 ***********************************************************************************************************/

int particleMetaStableIndex( Base const &a_particle ) {

    int metaStableIndex = 0;

    if( a_particle.isMetaStableAlias( ) ) {
        MetaStable const &object = (MetaStable const &) a_particle;
        metaStableIndex = object.metaStableIndex( );
    }

    return( metaStableIndex );
}

/* *********************************************************************************************************//**
 * Returns the meta-stable index if *a_index* is a PoPI::MetaStable alias; otherwise, 0 is returned. 
 *
 * @param a_pops                        [in]    The PoPs database to look up the particle.
 * @param a_index                       [in]    The index of the particle in *a_pops* whose meta-stable index value is returned.
 *
 * @return                                      The meta-stable index of the particle
 ***********************************************************************************************************/

int particleMetaStableIndex( Database const &a_pops, int a_index ) {

    Base const &object( a_pops.get<Base>( a_pops.final( a_index ) ) );

    return( particleMetaStableIndex( object ) );
}

/* *********************************************************************************************************//**
 * Returns the meta-stable index if *a_id* is a PoPI::MetaStable alias; otherwise, 0 is returned. 
 *
 * @param a_pops                        [in]    The PoPs database to look up the particle.
 * @param a_id                          [in]    The id of the particle in *a_pops* whose meta-stable index value is returned.
 *
 * @return                                      The meta-stable index of the particle
 ***********************************************************************************************************/

int particleMetaStableIndex( Database const &a_pops, std::string const &a_id ) {

    Base const &object( a_pops.get<Base>( a_pops.final( a_id ) ) );

    return( particleMetaStableIndex( object ) );
}

/* *********************************************************************************************************//**
 * A physical quantity can be a double, integer, fraction (e.g, '3/7') or a string. For all but string,
 * this functions returns a double representing the value of the physical quantity *a_physicalQuantity*.
 * For string, a throw is executed.
 *
 * @param a_physicalQuantity        [in]    The physical quantity whose value is returend.
 *
 * @return                                  A double value representing the physical quantity.
 ***********************************************************************************************************/

double getPhysicalQuantityAsDouble( PhysicalQuantity const &a_physicalQuantity ) {

    double value = 0.0;

    switch( a_physicalQuantity.Class( ) ) {
    case PQ_class::Double :
    case PQ_class::shell : {
        PQ_double const &pq_double = static_cast<PQ_double const &>( a_physicalQuantity );
        value = pq_double.value( ); }
        break;
    case PQ_class::integer : {
        PQ_integer const &pq_integer = static_cast<PQ_integer const &>( a_physicalQuantity );
        value = pq_integer.value( ); }
        break;
    default :
        throw Exception( "Cannot convert physical quantitiy to a double." );
    }

    return( value );
}

/* *********************************************************************************************************//**
 * If the suite *a_suite* as data, **getPhysicalQuantityAsDouble** is called on its first item; otherwise, 
 * a throw is executed. If *a_suite* is empty and *a_allowEmpty* is **true**, then *a_emptyValue* is returned.
 *
 * @param a_suite                   [in]    The suite whose first item's value is returned as a double.
 * @param a_allowEmpty              [in]    Determines act to follow when *a_suite* is empty.
 * @param a_emptyValue              [in]    The value to return if *a_suite* is empty and *a_allowEmpty* is **true**.
 *
 * @return                                  A double value representing the physical quantity.
 ***********************************************************************************************************/

double getPhysicalQuantityOfSuiteAsDouble( PQ_suite const &a_suite, bool a_allowEmpty, double a_emptyValue ) {

    if( a_suite.size( ) == 0 ) {
        if( a_allowEmpty ) return( a_emptyValue );
        throw Exception( "No physical quantitiy in Suite." );
    }

    return( getPhysicalQuantityAsDouble( *a_suite[0] ) );
}

/* *********************************************************************************************************//**
 * Breaks th components of a particles name into base, anti and quailier strings. Every id in GNDS PoPs can be of the
 * form base["_anti"][qualifier] where both "_anti" and qualifier are optional. For example, an electron is represented
 * as "e-" and an anti-electron (i.e., positron) "e-_anti". Qualifiers are endings imbedded by "{" and "}". For example,
 * "H{1s1/2}" has the base "H" with quailifier "1s1/2" where the "{" and "}" have been stripped from the quailifier.
 *
 * @param a_id                      [in]    The base id for *a_id*.
 * @param a_anti                    [in]    A std::string to be filled with "_anti" if particle is an anti-particle and an "" otherwise.
 * @param a_qualifier               [in]    A pointer to a std::string that will be filled with the qualifer characters.
 *
 * @return                                  The base id for the particle.
 ***********************************************************************************************************/

std::string baseAntiQualifierFromID( std::string const &a_id, std::string &a_anti, std::string *a_qualifier) {

    std::size_t curlyBraketPosition = a_id.find( "{" );
    std::string base = a_id.substr( 0, curlyBraketPosition );

    a_anti = "";
    if( a_qualifier != nullptr ) *a_qualifier = "";

    if( curlyBraketPosition != std::string::npos ) {
        if( a_id.back( ) != '}' ) throw Exception( "Invalid quaifier string in id '" + a_id + "'." );
        base = a_id.substr( 0, curlyBraketPosition );
        if( a_qualifier != nullptr ) {
            *a_qualifier = a_id.substr( curlyBraketPosition + 1, a_id.size( ) - curlyBraketPosition - 2 );
        } }
    else if( a_id.find( "}" ) != std::string::npos ) {
        throw Exception( "Invalid quaifier string in id '" + a_id + "'." );
    }

    std::size_t anti_position = base.find( IDs::anti );
    if( anti_position != std::string::npos ) {
        a_anti = base.substr( anti_position );
        base = base.substr( 0, anti_position );
        if( a_anti != IDs::anti ) throw Exception( "Invalid anti string in id '" + a_id + "'." );
    }

    return( base );
}

/*! \class Exception
 * Exception class for all PoPI exceptions thrown by PoPI functions.
 */

/* *********************************************************************************************************//**
 * @param a_message         [in]     The message that the function what() will return.
 ***********************************************************************************************************/

Exception::Exception( std::string const & a_message ) :
        std::runtime_error( a_message ) {

}

}
