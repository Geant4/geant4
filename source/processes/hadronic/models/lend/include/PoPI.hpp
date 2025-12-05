/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef PoPI_hpp_included
#define PoPI_hpp_included 1

#include <string>
#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <stdexcept>
#include <typeinfo>
#include <fstream>
#include <exception>
#include <utility>
#include <stddef.h>

#include <LUPI.hpp>
#include <HAPI.hpp>

namespace PoPI {

#define PoPI_AMU2MeV_c2 931.494028
#define PoPI_electronMass_MeV_c2 0.5109989461

#define PoPI_formatVersion_0_1_Chars "0.1"
#define PoPI_formatVersion_1_10_Chars "1.10"
#define PoPI_formatVersion_2_0_LLNL_3_Chars "2.0.LLNL_3"

#define PoPI_PoPsChars "PoPs"

#define PoPI_idChars "id"
#define PoPI_symbolChars "symbol"
#define PoPI_chemicalElementsChars "chemicalElements"
#define PoPI_chemicalElementChars "chemicalElement"
#define PoPI_isotopesChars "isotopes"
#define PoPI_isotopeChars "isotope"
#define PoPI_gaugeBosonChars "gaugeBoson"
#define PoPI_leptonChars "lepton"
#define PoPI_baryonChars "baryon"
#define PoPI_nuclidesChars "nuclides"
#define PoPI_nuclideChars "nuclide"
#define PoPI_nucleusChars "nucleus"
#define PoPI_unorthodoxChars "unorthodox"
#define PoPI_aliasesChars "aliases"

/*! \enum Particle_class
 * This enum represents the various type of allowed particle types.
 */

enum class Particle_class { nuclide,                /**< Specifies that the particle is a nuclide. */
                            nucleus,                /**< Specifies that the particle is a nucleus. */
                            gaugeBoson,             /**< Specifies that the particle is a gauge boson. */
                            lepton,                 /**< Specifies that the particle is a lepton. */
                            baryon,                 /**< Specifies that the particle is a baryon. */
                            nuclideMetaStable,      /**< Specifies that the particle is a nuclide meta-stable alias. */
                            nucleusMetaStable,      /**< Specifies that the particle is a nucleus meta-stable alias. */
                            TNSL,                   /**< Specifies that the particle is a TNSL target. Currently not used. */
                            ENDL_fissionProduct,    /**< Specifies that the particle is an ENDL fissiont product (e.g., 99120, 99125). */
                            unorthodox,             /**< Specifies that the particle is an unorthodox. */
                            alias,                  /**< Specifies that the particle is a alias. */
                            chemicalElement,        /**< Specifies that the particle is a chemicalElement. */
                            isotope,                /**< Specifies that the particle is a isotope. */
                            unknown                 /**< Specifies that the particle is a unknown. */ };

#define PoPI_massChars "mass"
#define PoPI_spinChars "spin"
#define PoPI_parityChars "parity"
#define PoPI_chargeChars "charge"
#define PoPI_halflifeChars "halflife"

#define PoPI_doubleChars "double"
#define PoPI_integerChars "integer"
#define PoPI_fractionChars "fraction"
#define PoPI_stringChars "string"
#define PoPI_shellChars "shell"
#define PoPI_decayDataChars "decayData"
#define PoPI_gammaDecayDataChars "gammaDecayData"

#define PoPI_decayModeElectroMagnetic "electroMagnetic"

#define PoPI_formatChars "format"
#define PoPI_labelChars "label"
#define PoPI_indexChars "index"
#define PoPI_pidChars "pid"
#define PoPI_nameChars "name"
#define PoPI_versionChars "version"
#define PoPI_aliasChars "alias"
#define PoPI_metaStableChars "metaStable"
#define PoPI_particleChars "particle"
#define PoPI_discreteChars "discrete"
#define PoPI_continuumChars "continuum"

/*! \enum PQ_class
 * This enum represents the various type of allowed physcial quantity types.
 */

enum class PQ_class { Double,       /**< Specifies that the physcial quantity is a double. */
                      integer,      /**< Specifies that the physcial quantity is an integer. */
                      fraction,     /**< Specifies that the physcial quantity is a fraction. */
                      string,       /**< Specifies that the physcial quantity is a string. */
                      shell         /**< Specifies that the physcial quantity is a shell. */ };

/*! \enum SpecialParticleID_mode
 * This enum specifies how the light charged particle ids are handled. The light charged particles ids are familiarly known as
 * p, d, t, h and a.
 */

enum class SpecialParticleID_mode { familiar,   /**< Treat ids as the familiar p, d, t, h and a. */
                                    nuclide,    /**< Treat ids as the familiar p, d, t, h and a as h1, h2, h3, he3 and he4, respectively. */
                                    nucleus     /**< Treat ids as the familiar p, d, t, h and a as H1, H2, H3, He3 and He4, respectively. */ };

class NuclideGammaBranchStateInfos;
class Base;
class SymbolBase;
class Decay;
class DecayMode;
class DecayData;
class Particle;
class MetaStable;
class Alias;
class Baryon;
class GaugeBoson;
class Lepton;
class Nuclide;
class Nucleus;
class Unorthodox;

class Isotope;
class ChemicalElement;
class Database;

void appendXMLEnd( std::vector<std::string> &a_XMLList, std::string const &a_label );

int particleZ( Base const &a_particle, bool a_isNeutronProtonANucleon = false );
int particleZ( Database const &a_pops, std::size_t a_index, bool a_isNeutronProtonANucleon = false );
int particleZ( Database const &a_pops, std::string const &a_id, bool a_isNeutronProtonANucleon = false );

int particleA( Base const &a_particle, bool a_isNeutronProtonANucleon = false );
int particleA( Database const &a_pops, std::size_t a_index, bool a_isNeutronProtonANucleon = false );
int particleA( Database const &a_pops, std::string const &a_id, bool a_isNeutronProtonANucleon = false );

int particleZA( Base const &a_particle, bool a_isNeutronProtonANucleon = false );
int particleZA( Database const &a_pops, std::size_t a_index, bool a_isNeutronProtonANucleon = false );
int particleZA( Database const &a_pops, std::string const &a_id, bool a_isNeutronProtonANucleon = false );

int particleMetaStableIndex( Base const &a_particle );
int particleMetaStableIndex( Database const &a_pops, std::size_t a_index );
int particleMetaStableIndex( Database const &a_pops, std::string const &a_id );

std::string specialParticleID( SpecialParticleID_mode a_mode, std::string const &a_id );
bool compareSpecialParticleIDs( std::string const &a_id1, std::string const &a_id2 );

struct IDs {
    static std::string const photon;
    static std::string const electron;
    static std::string const neutron;
    static std::string const proton;
    static std::string const familiarPhoton;
    static std::string const familiarDeuteron;
    static std::string const familiarTriton;
    static std::string const familiarHelion;
    static std::string const familiarAlpha;
    static std::string const FissionProductENDL99120;
    static std::string const FissionProductENDL99125;
    static std::string const anti;
};

struct Intids {
    static int constexpr neutron = 1020000000;
    static int constexpr photon = 1000000000;
    static int constexpr electron = 1010000000;
    static int constexpr FissionProductENDL99120 = 1990099120;
    static int constexpr FissionProductENDL99125 = 1990099125;
};

extern std::map<std::string, std::string> supportedNucleusAliases;

typedef std::vector<Base *> ParticleList;
typedef std::vector<SymbolBase *> SymbolList;

/*
============================================================
======================== Exception =========================
============================================================
*/

class Exception : public std::runtime_error {

    public :
        explicit Exception( std::string const &a_message );

};

/*
============================================================
====================== ParseIntidInfo ======================
============================================================
*/

class ParseIntidInfo {

    private:
        int m_intid;                        /**< The intid for the rest of the data. */
        Particle_class m_family;            /**< The family of the particle. */
        bool m_isAnti;                      /**< **true** if particle is an anti-particle and **false** otherwise. */
// The following are for nuclear like particles.
        bool m_isNuclear;                   /**< *true* if the particle is a nuclear particle and *false* otherwise. */
        int m_AAA;                          /**< For a nuclear particle, its AAA value (i.e., atomic mass number). */
        int m_ZZZ;                          /**< For a nuclear particle, its ZZZ value (i.e., atomic number). */
        int m_III;                          /**< For a nuclear particle, its III value (nuclear excitation level index or meta-stable index. */
        int m_nuclearLevelIndex;            /**< For a nuclear particle, nuclear excitation level index. */
        int m_metaStableIndex;              /**< For a nuclear meta-stable particle, its meta-stable index. */
// The following are for leptons.
        int m_generation;                   /**< For a lepton, its generation. */
        bool m_isNeutrino;                  /**< For a lepton, **true** if leption is a neutrino and **false** otherwise. */
// The follow are for baryons.
        int m_baryonGroup;                  /**< For a baryon, its baryon group. */
        int m_baryonId;                     /**< For a baryon, its id within a baryon group. */
// Other data.
        int m_familyId;                     /**< For non-nuclear particles, the particle's indentifier within its family. */

    public:
        ParseIntidInfo( int a_intid, bool a_GRIN_mode = false );

        int intid( ) { return( m_intid ); }                     /**< Returns the value of the *m_intid* member. */
        Particle_class family( ) { return( m_family ); }        /**< Returns the value of the *m_family* member. */
        bool isAnti( ) { return( m_isAnti ); }                  /**< Returns the value of the *m_isAnti * member. */

        bool isNuclear( ) { return( m_isNuclear ); }            /**< Returns the value of the *m_isNuclear* member. */
        int AAA( ) { return( m_AAA ); }                         /**< Returns the value of the *m_AAA* member. */
        int ZZZ( ) { return( m_ZZZ ); }                         /**< Returns the value of the *m_ZZZ* member. */
        int III( ) { return( m_III ); }                         /**< Returns the value of the *m_III* member. */
        bool isNuclearMetaStable( ) { return( ( m_family == Particle_class::nuclideMetaStable ) || ( m_family == Particle_class::nucleusMetaStable ) ); }
                                                                /**< Returns **true** if particle is a nuclear meta-stable alias. */
        int metaStableIndex( ) { return( m_metaStableIndex ); } /**< Returns the value of the *m_metaStableIndex* member. */

        int generation( ) { return( m_generation ); }           /**< Returns the value of the *m_generation* member. */
        bool isNeutrino( ) { return( m_isNeutrino ); }          /**< Returns the value of the *m_isNeutrino* member. */

        int baryonGroup( ) { return( m_baryonGroup ); }         /**< Returns the value of the *m_baryonGroup* member. */
        int baryonId( ) { return( m_baryonId ); }               /**< Returns the value of the *m_baryonId * member. */

        int familyId( ) { return( m_familyId ); }               /**< Returns the value of the *m_familyId* member. */

        std::string id( );
};

/*
============================================================
======================= ParseIdInfo ========================
============================================================
*/

class ParseIdInfo{

    private:
        bool m_isSupported;                                 /**< If **true** the particle's id was parsed and the other member of *this* are valid. Otherwise, parsing of the particle's id is currently not supported. */
        std::string m_id;                                   /**< The id for the particles. */
        bool m_isNuclear;                                   /**< **true** if particle is a valid nuclear name (i.e., nuclide, nucleus of meta-stable) and **false** otherwise. */
        bool m_isNucleus;                                   /**< **true** if particle is a nucleus and **false** otherwise. */
        bool m_isChemicalElement;                           /**< **true** if id is only a chemical element symbol. */
        bool m_isAnti;                                      /**< **true** if particle is an anti-particle and **false** otherwise. */
        bool m_isMetaStable;                                /**< **true** if particle is a meta-stable alias and **false** otherwise. */
        std::string m_symbol;                               /**< The chemical element symbol part of *m_id*. This will always be the nuclide symbol even if *m_id* is for a nucleus. */
        int m_Z;                                            /**< The atomic number of the *m_id*. */
        int m_A;                                            /**< The atomic mass number of the *m_id*. */
        int m_index;                                        /**< The nuclear level index of *m_id*. */
        std::string m_qualifier;

        std::string boolToString( bool a_value, std::string const &a_prefix ) const;

    public:
        ParseIdInfo( std::string const &a_id );

        bool isSupported( ) { return( m_isSupported ); }                /**< Returns the value of the *m_isSupported. */
        std::string const &Id( ) { return( m_id ); }                    /**< Returns a reference to the *m_id* member. */
        bool isNuclear( ) { return( m_isNuclear ); }                    /**< Returns the value of the *m_isNuclear* member. */
        bool isNucleus( ) { return( m_isNucleus ); }                    /**< Returns the value of the *m_isNucleus* member. */
        bool isChemicalElement( ) { return( m_isChemicalElement ); }    /**< Returns the value of the *m_isChemicalElement* member. */
        bool isAnti( ) { return( m_isAnti ); }                          /**< Returns the value of the *m_isAnti* member. */
        std::string const &symbol( ) { return( m_symbol); }             /**< Returns a reference to the *m_symbol* member. */
        int Z( ) { return( m_Z ); }                                     /**< Returns the value of the *m_Z* member. */
        int A( ) { return( m_A ); }                                     /**< Returns the value of the *m_A* member. */
        int index( ) { return( m_index ); }                             /**< Returns the value of the *m_index* member. */
        std::string const &qualifier( ) { return( m_qualifier ); }      /**< Returns a reference to the *m_qualifier* member. */

        void print( bool a_terse, std::string const &a_indent = "" ) const ;
};

/*! \class Suite
 * This is the base class for all suite like members.
 */

/*
============================================================
========================== Suite ===========================
============================================================
*/

template <class T, class T2>
class Suite {

    private:
        std::string m_moniker;                                  /**< The moniker (i.e., name) of the suite. */
        std::vector<T *> m_items;                               /**< The list of all items in the suite. */

    public:
        Suite( std::string const &a_moniker ) : m_moniker( a_moniker ) { }
        ~Suite( );
        void appendFromParentNode( HAPI::Node const &a_node, Database *a_DB, T2 *a_parent );
        void appendFromParentNode2( HAPI::Node const &a_node, T2 *a_parent );

        std::string::size_type size( void ) const { return( m_items.size( ) ); }        /**< Returns the number of items in the suite. */
        T &operator[]( std::size_t a_index ) const { return( *m_items[a_index] ); }     /**< Returns the item at index *a_index*. */
        std::string const &moniker( void ) { return( m_moniker ); }                     /**< Returns the value of the *m_moniker* member. */

        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

template <class T, class T2>
Suite<T, T2>::~Suite( ) {

    std::string::size_type i1, _size = m_items.size( );

    for( i1 = 0; i1 < _size; ++i1 ) delete m_items[i1];

// Ask Adam why next line does not work.
//    for( std::vector<T *>::iterator iter = m_items.begin( ); iter != m_items.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Adds the children of *a_node* to the suite and to *a_DB*.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_DB              [in]    The **PoPI::Database** instance to add the constructed items to.
 * @param a_parent          [in]    The parent suite that will contain *this*.
 ***********************************************************************************************************/

template <class T, class T2>
void Suite<T, T2>::appendFromParentNode( HAPI::Node const &a_node, Database *a_DB, T2 *a_parent ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        T *item = new T( child, a_DB, a_parent );
        m_items.push_back( item );
    }
}

/* *********************************************************************************************************//**
 * Adds the children of *a_node* to the suite.
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_parent          [in]    The parent suite that will contain *this*.
 ***********************************************************************************************************/

template <class T, class T2>
void Suite<T, T2>::appendFromParentNode2( HAPI::Node const &a_node, T2 *a_parent ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        T *item = new T( child, a_parent );
        m_items.push_back( item );
    }
}

/* *********************************************************************************************************//**
 * Creates an XML representation of the suite.
 *
 * @param a_XMLList         [in]    The list the XML lines are added to.
 * @param a_indent1         [in]    The amount to indent the XML text.
 ***********************************************************************************************************/

template <class T, class T2>
void Suite<T, T2>::toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const {

    std::string::size_type _size = m_items.size( );
    std::string indent2 = a_indent1 + "  ";

    if( _size == 0 ) return;

    std::string header = a_indent1 + "<" + m_moniker + ">";
    a_XMLList.push_back( std::move( header ) );
    for( std::string::size_type i1 = 0; i1 < _size; ++i1 ) m_items[i1]->toXMLList( a_XMLList, indent2 );

    appendXMLEnd( a_XMLList, m_moniker );
}

/*
============================================================
===================== PhysicalQuantity =====================
============================================================
*/

class PhysicalQuantity {

    private:
        PQ_class m_class;                           /**< The class for the physical quanity. */
        std::string m_tag;                          /**< The name of the physical quanity. */
        std::string m_label;                        /**< The label for the physical quanity. */
        std::string m_valueString;                  /**< The string value of the physical quanity. */
        std::string m_unit;                         /**< The unit of the physical quanity. */

    public:
        PhysicalQuantity( HAPI::Node const &a_node, PQ_class a_class );
        virtual ~PhysicalQuantity( );

        PQ_class Class( void ) const { return( m_class ); }                         /**< Returns the value of the *m_class* member. */
        std::string const &tag( void ) const { return( m_tag ); }                   /**< Returns the value of the *m_tag* member. */
        std::string const &label( void ) const { return( m_label ); }               /**< Returns the value of the *m_label* member. */
        std::string const &valueString( void ) const { return( m_valueString ); }   /**< Returns the value of the *valueString* member. */
        std::string const &unit( void ) const { return( m_unit ); }                 /**< Returns the value of the *m_unit* member. */

        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
        virtual std::string valueToString( void ) const = 0;
};

/*
============================================================
========================= PQ_double ========================
============================================================
*/

class PQ_double : public PhysicalQuantity {

    private:
        double m_value;                         /**< The double value of the physical quanity. */
        void initialize( );

    public:
        PQ_double( HAPI::Node const &a_node );
        PQ_double( HAPI::Node const &a_node, PQ_class a_class );
        virtual ~PQ_double( );

        double value( void ) const { return( m_value ); }                           /**< Returns the value of the *m_value* member. */
        double value( char const *a_unit ) const ;
        double value( std::string const &a_unit ) const { return( value( a_unit.c_str( ) ) ); }
                                                                                    /**< Returns the value of the *m_value* member in units of *a_unit*. */
        virtual std::string valueToString( void ) const ;
};

/*
============================================================
========================= PQ_integer =======================
============================================================
*/

class PQ_integer : public PhysicalQuantity {

    private:
        int m_value;                            /**< The integer value of the physical quanity. */

    public:
        PQ_integer ( HAPI::Node const &a_node );
        virtual ~PQ_integer( );

        int value( void ) const { return( m_value ); }                              /**< Returns the value of the *m_value* member. */
        int value( char const *a_unit ) const ;
        int value( std::string const &a_unit ) const { return( value( a_unit.c_str( ) ) ); }
        virtual std::string valueToString( void ) const ;
};

/*
==============================================================
========================= PQ_fraction ========================
==============================================================
*/

class PQ_fraction : public PhysicalQuantity {

    public:
        PQ_fraction( HAPI::Node const &a_node );
        virtual ~PQ_fraction( );

        std::string value( void ) const ;
        std::string value( char const *a_unit ) const ;
        std::string value( std::string const &a_unit ) const { return( value( a_unit.c_str( ) ) ); }
                                                                                    /**< Returns the value of the *m_value* member in units of *a_unit*. */
        virtual std::string valueToString( void ) const ;
};

/*
============================================================
========================= PQ_string ========================
============================================================
*/

class PQ_string : public PhysicalQuantity {

    public:
        PQ_string( HAPI::Node const &a_node );
        virtual ~PQ_string( );

        std::string value( void ) const { return( valueString( ) ); }               /**< Returns the value returned by calling the **valueString** methods. */
        std::string value( char const *a_unit ) const ;
        std::string value( std::string const &a_unit ) const { return( value( a_unit.c_str( ) ) ); }
                                                                                    /**< Returns the value of the *m_value* member in units of *a_unit*. */
        virtual std::string valueToString( void ) const ;
};

/*
============================================================
========================= PQ_shell =========================
============================================================
*/

class PQ_shell : public PQ_double {

    public:
        PQ_shell( HAPI::Node const &a_node );
        ~PQ_shell( );
};

/*
============================================================
========================= PQ_suite =========================
============================================================
*/

class PQ_suite : public std::vector<PhysicalQuantity *> {

    private:
        std::string m_label;

    public:
        PQ_suite( HAPI::Node const &a_node );
        ~PQ_suite( );

        std::string &label( void ) { return( m_label ); }

        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
================== NuclideGammaBranchInfo ==================
============================================================
*/

class NuclideGammaBranchInfo {

    private:
        double m_probability;                               /**< The probability that the level decays to state *m_residualState*. */
        double m_photonEmissionProbability;                 /**< The conditional probability the the decay emitted a photon. */
        double m_gammaEnergy;                               /**< The energy of the emitted photon. */
        std::string m_residualState;                        /**< The state the residual is left in after photon decay. */

    public:
        NuclideGammaBranchInfo( double a_probability, double a_photonEmissionProbability, double a_gammaEnergy, std::string const &a_residualState );
        NuclideGammaBranchInfo( NuclideGammaBranchInfo const &a_nuclideGammaBranchInfo );

        double probability( ) const { return( m_probability ); }                        /**< Returns the value of the *m_probability* member. */
        double photonEmissionProbability( ) const { return( m_photonEmissionProbability ); }
                                                                                        /**< Returns the value of the *m_photonEmissionProbability* member. */
        double gammaEnergy( ) const { return( m_gammaEnergy ); }                        /**< Returns the value of the *m_gammaEnergy* member. */
        std::string const &residualState( ) const { return( m_residualState ); }        /**< Returns the value of the *m_residualState* member. */
};      

/*
==============================================================
================= NuclideGammaBranchStateInfo ================
==============================================================
*/

class NuclideGammaBranchStateInfo {

    private:
        std::string m_state;                                    /**< The inital state the decay starts from. */
        int m_intid;                                            /**< The intid for the inital state. */
        std::string m_kind;                                     /**< The kind of the particle. Currently can be 'discrete' or 'continuum'. */
        double m_nuclearLevelEnergy;                            /**< The nuclear level excitation energy of the level (state). */
        double m_nuclearLevelEnergyWidth;                       /**< This is 0.0 except for GRIN realized continuum levels where this is the energy width from this level to the next higher level. */
        bool m_derivedCalculated;                               /**< For internal use to determine if other members have been set or not. */
        double m_multiplicity;                                  /**< The average number of photons emitted when transitioning from the initial to the final state. Data derived from m_branches data. */
        double m_averageGammaEnergy;                            /**< The average energy per decay from the initial to the final state. Data derived from m_branches data. */
        std::vector<NuclideGammaBranchInfo> m_branches;

    public:
        NuclideGammaBranchStateInfo( std::string const &a_state, int a_intid, std::string const &a_kind, double a_nuclearLevelEnergy );

        std::string const &state( ) const { return( m_state ); }                /**< Returns the value of the *m_state* member. */
        int intid( ) const { return( m_intid ); }                               /**< Returns the value of the *m_intid* member. */
        std::string const &kind( ) const { return( m_kind ); }                  /**< Returns the value of the *m_kind* member. */
        double nuclearLevelEnergy( ) const { return( m_nuclearLevelEnergy ); }  /**< Returns the value of the *m_nuclearLevelEnergy* member. */
        double nuclearLevelEnergyWidth( ) const { return( m_nuclearLevelEnergyWidth ); }
                                                                                /**< Returns the value of the *m_nuclearLevelEnergyWidth* member. */
        void setNuclearLevelEnergyWidth( double a_nuclearLevelEnergyWidth ) { m_nuclearLevelEnergyWidth = a_nuclearLevelEnergyWidth; }
                                                                                /**< Set the value of the *m_nuclearLevelEnergyWidth* member to *a_nuclearLevelEnergyWidth*. */
        bool derivedCalculated( ) const { return( m_derivedCalculated ); }      /**< Returns the value of the *m_derivedCalculated* member. */
        double multiplicity( ) const { return( m_multiplicity ); }              /**< Returns the value of the *m_multiplicity* member. */
        double averageGammaEnergy( ) const { return( m_averageGammaEnergy ); }  /**< Returns the value of the *m_averageGammaEnergy* member. */
        std::vector<NuclideGammaBranchInfo> const &branches( ) const { return( m_branches ); }
                                                                                /**< Returns a reference to the *m_branches* member. */

        void add( NuclideGammaBranchInfo const &a_nuclideGammaBranchInfo );
        void calculateDerivedData( NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos );
};

/*
==============================================================
================ NuclideGammaBranchStateInfos ================
==============================================================
*/

class NuclideGammaBranchStateInfos {

    private:
        std::vector<NuclideGammaBranchStateInfo *> m_nuclideGammaBranchStateInfos;

    public:
        NuclideGammaBranchStateInfos( );
        ~NuclideGammaBranchStateInfos( );

        std::size_t size( ) const { return( m_nuclideGammaBranchStateInfos.size( ) ); }
        NuclideGammaBranchStateInfo *operator[]( std::size_t a_index ) { return( m_nuclideGammaBranchStateInfos[a_index] ); }
        NuclideGammaBranchStateInfo const *operator[]( std::size_t a_index ) const { return( m_nuclideGammaBranchStateInfos[a_index] ); }
        std::vector<NuclideGammaBranchStateInfo *> &nuclideGammaBranchStateInfos( ) { return( m_nuclideGammaBranchStateInfos ); }
        void add( NuclideGammaBranchStateInfo *a_nuclideGammaBranchStateInfo );
        NuclideGammaBranchStateInfo *find( std::string const &a_state );
        NuclideGammaBranchStateInfo const *find( std::string const &a_state ) const ;
};

/*
============================================================
=========================== Base ===========================
============================================================
*/

class Base {

    private:
        std::string m_id;                               /**< The **PoPs** id for the particle or **PoPs** symbol for a chemicalElement or isotope. */
        Particle_class m_class;                         /**< The **Particle_class** for the particle, chemicalElement or isotope. */
        std::size_t m_index;                            /**< The for the particle, chemicalElement or isotope. */
        int m_intid;                                    /**< The unique integer id for a particle or a meta-stable alias. For a non meta-stable alias, an isotope or chemical element, this is -1. */

        void setIntid( int a_intid ) { m_intid = a_intid; }                                 /**< Sets the value of the *m_intid* member to *a_intid*. */

    public:
        Base( std::string const &a_id, Particle_class a_class );
        Base( HAPI::Node const &a_node, std::string const &a_label, Particle_class a_class );
        virtual ~Base( );

        std::string const &ID( void ) const { return( m_id ); }                             /**< Returns a *const* reference to the *m_id* member of *this*. */
        std::size_t index( void ) const { return( m_index ); }                              /**< Returns the value of the *m_index* member of *this*. */
        void setIndex( std::size_t a_index ) { m_index = a_index; }                         /**< Sets the value of the *m_index* member of *this* to *a_index*. */
        int intid( ) const { return( m_intid ); }                                           /**< Returns the value of the *m_intid* member. */
        Particle_class Class( void ) const { return( m_class ); }                           /**< Returns the value of the *m_class* member of *this*. */
        virtual bool isParticle( ) const { return( true ); }                                /**< Returns **true** if *this* is a **Particle** and **false** it *this* is a **ChemicalElement** or **Isotope** instance. */
        bool isAlias( void ) const { return( ( m_class == Particle_class::alias ) || isMetaStableAlias( ) ); }
                                                                                            /**< Returns **true** if *this* is an **Alias** or **MetaStable** instance and **false** otherwise. */
        bool isMetaStableAlias( void ) const { return( ( m_class == Particle_class::nuclideMetaStable ) || ( m_class ==  Particle_class::nucleusMetaStable ) ); }
                                                                                            /**< Returns **true** if *this* is a **MetaStable** instance and **false** otherwise. */

        bool isGaugeBoson( ) const { return( m_class == Particle_class::gaugeBoson ); }     /**< Returns **true** if *this* is a **GaugeBoson** instance and **false** otherwise. */
        bool isLepton( ) const { return( m_class == Particle_class::lepton ); }             /**< Returns **true** if *this* is a **Lepton** instance and **false** otherwise. */
        bool isBaryon( ) const { return( m_class == Particle_class::baryon ); }             /**< Returns **true** if *this* is a **Baryon** instance and **false** otherwise. */
        bool isUnorthodox( ) const { return( m_class == Particle_class::unorthodox ); }     /**< Returns **true** if *this* is a **Unorthodox** instance and **false** otherwise. */
        bool isNucleus( ) const { return( m_class == Particle_class::nucleus ); }           /**< Returns **true** if *this* is a **Nucleus** instance and **false** otherwise. */
        bool isNuclide( ) const { return( m_class == Particle_class::nuclide ); }           /**< Returns **true** if *this* is a **Nuclide** instance and **false** otherwise. */
        bool isIsotope( ) const { return( m_class == Particle_class::isotope ); }           /**< Returns **true** if *this* is a **Isotope** instance and **false** otherwise. */
        bool isChemicalElement( ) const { return( m_class == Particle_class::chemicalElement ); }
                                                                                            /**< Returns **true** if *this* is a **ChemicalElement** instance and **false** otherwise. */

        friend MetaStable;
        friend Alias;
        friend Baryon;
        friend GaugeBoson;
        friend Lepton;
        friend Nucleus;
        friend Nuclide;
        friend Unorthodox;
};

/*
============================================================
========================== IDBase ==========================
============================================================
*/

class IDBase : public Base {

    public:
        IDBase( std::string const &a_id, Particle_class a_class );
        IDBase( HAPI::Node const &a_node, Particle_class a_class );
        virtual ~IDBase( );       // BRB This should be virtual but I cannot get it to work without crashing.

        std::size_t addToDatabase( Database *a_DB );
        double massValue2( Database const &a_DB, std::string const &a_unit ) const ;
};

/*
============================================================
======================== SymbolBase ========================
============================================================
*/

class SymbolBase : public Base {

    public:
        SymbolBase( HAPI::Node const &a_node, Particle_class a_class );
        ~SymbolBase( );

        std::string const &symbol( ) const { return( ID( ) ); }                             /**< Returns the value of the symbol. */

        std::size_t addToSymbols( Database *a_DB );
        bool isParticle( ) const { return( false ); }
};

/*
============================================================
========================= Product ==========================
============================================================
*/

class Product {

    private:
        int m_id;
        std::string m_pid;
        std::string m_label;

    public:
        Product( HAPI::Node const &a_node, Decay *a_DB );
        ~Product( );

        int ID( ) const { return( m_id ); }
        std::string const &pid( ) const { return( m_pid ); }
        std::string const &label( ) const { return( m_label ); }

        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
========================== Decay ===========================
============================================================
*/

class Decay {

    private:
        int m_index;
        std::string m_mode;
        bool m_complete;
        Suite<Product, Decay> m_products;

    public:
        Decay( HAPI::Node const &a_node, DecayMode const *a_decayMode );
        ~Decay( );

        int index( void ) const { return( m_index ); }
        std::string const &mode( ) const { return( m_mode ); }
        bool complete( ) const { return( m_complete ); }
        Suite<Product, Decay> const &products( void ) const { return( m_products ); }
        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
======================== DecayMode =========================
============================================================
*/

class DecayMode {

    private:
        std::string m_label;
        std::string m_mode;
        PQ_suite m_probability;
        PQ_suite m_photonEmissionProbabilities;
        Suite<Decay, DecayMode> m_decayPath;

    public:
        DecayMode( HAPI::Node const &a_node, DecayData const *a_decayData );
        ~DecayMode( );

        std::string const &label( ) const { return( m_label ); }
        std::string const &mode( ) const { return( m_mode ); }
        PQ_suite const &probability( ) const { return( m_probability ); }
        PQ_suite const &photonEmissionProbabilities( ) const { return( m_photonEmissionProbabilities ); }
        Suite<Decay, DecayMode> const &decayPath( ) const { return( m_decayPath ); }

        void calculateNuclideGammaBranchStateInfo( PoPI::Database const &a_pops, NuclideGammaBranchStateInfo &a_nuclideGammaBranchStateInfo ) const ;
        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
======================== DecayData =========================
============================================================
*/

class DecayData {

    private:
        Suite<DecayMode, DecayData> m_decayModes;

    public:
        DecayData( HAPI::Node const &a_node );
        ~DecayData( );

        Suite<DecayMode, DecayData> const &decayModes( void ) const { return( m_decayModes ); }

        void calculateNuclideGammaBranchStateInfo( PoPI::Database const &a_pops, NuclideGammaBranchStateInfo &a_nuclideGammaBranchStateInfo ) const ;
        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
====================== GammaDecayData ======================
============================================================
*/

class GammaDecayData {

    private:
        std::string m_kind;                                     /**< The kind of the particle. Currently can be 'discrete' or 'continuum' but may be 'experimental', 'evaluated' or 'modelled' in the future. */
        int m_rows;                                             /**< The number of nuclides listed. */
        int m_columns;                                          /**< The number of items in a row. */
        std::vector<std::string> m_ids;                         /**< The list of nulcides. */
        std::vector<double> m_probabilities;                    /**< The list of probabilities for each nuclide. This must sum to 1. */
        std::vector<double> m_photonEmissionProbabilities;      /**< The list of photon emission probabilities for each nuclide. */

    public:
        GammaDecayData( HAPI::Node const &a_node );
        ~GammaDecayData( );

        std::string const &kind( ) const { return( m_kind ); }              /**< Returns a const reference to the *m_kind* member. */
        int rows( ) const { return( m_rows ); }
        int colunms( ) const { return( m_columns ); }
        std::vector<std::string> const &ids( ) const { return( m_ids ); }
        std::vector<double> const &probabilities( ) const { return( m_probabilities ); }
        std::vector<double> const &photonEmissionProbabilities( ) const { return( m_photonEmissionProbabilities ); }

        void calculateNuclideGammaBranchStateInfo( PoPI::Database const &a_pops, NuclideGammaBranchStateInfo &nuclideGammaBranchStateInfo ) const ;
};

/*
============================================================
========================= Particle =========================
============================================================
*/

class Particle : public IDBase {

    private:
        std::string m_baseId;                           /**< The base part of the id (i.e., without the anti and quailifier). */
        std::string m_family;                           /**< The family of the particle. */
        std::string m_anti;                             /**< The string "_anti" if particle is an anti-particle and an empty string otherwise. */
        int m_hasNucleus;                               /**< Indicates if the particle is or contains a nucleus. 0 = no, -1 = yes and 1 = is nucleus. */
        PQ_suite m_mass;                                /**< A suite storing the mass physical quantities for the particle. */
        PQ_suite m_spin;                                /**< A suite storing the spin physical quantities for the particle. */
        PQ_suite m_parity;                              /**< A suite storing the parity physical quantities for the particle. */
        PQ_suite m_charge;                              /**< A suite storing the charge physical quantities for the particle. */
        PQ_suite m_halflife;                            /**< A suite storing the halflife physical quantities for the particle. */
        DecayData m_decayData;                          /**< Stores the decay data for the particle. */

        void setHasNucleus( bool a_hasNucleus ) { m_hasNucleus = a_hasNucleus; }

    public:
        Particle( HAPI::Node const &a_node, Particle_class a_class, std::string const &a_family, int a_hasNucleus = 0 );
        virtual ~Particle( );

        std::string const &baseId( void ) const { return( m_baseId ); }         /**< Returns a *const* reference to the *m_baseId* member. */
        std::string const &family( void ) const { return( m_family ); }         /**< Returns a *const* reference to the *m_family* member. */
        bool isAnti( ) const { return( m_anti == IDs::anti ); }                              /**< Returns the value of the *m_anti* member. */
        int hasNucleus( void ) const { return( m_hasNucleus ); }                /**< Returns the value of the *m_hasNucleus* member. */

        virtual PQ_suite const &mass( void ) const { return( m_mass ); }        /**< Returns a *const* reference to the *m_mass* member. */
        virtual double massValue( char const *a_unit ) const ;
        double massValue( std::string const &a_unit ) const { return( massValue( a_unit.c_str( ) ) ); }
                                                                                /**< Returns the value of massValue( a_unit.c_str( ) ). */

        PQ_suite const &spin( ) const { return( m_spin ); }                     /**< Returns a *const* reference to the *m_spin* member. */
        PQ_suite const &parity( ) const { return( m_parity ); }                 /**< Returns a *const* reference to the *m_parity* member. */
        PQ_suite const &charge( ) const { return( m_charge ); }                 /**< Returns a *const* reference to the *m_charge* member. */
        PQ_suite const &halflife( ) const { return( m_halflife ); }             /**< Returns a *const* reference to the *m_halflife* member. */
        DecayData const &decayData( ) const { return( m_decayData ); }          /**< Returns a *const* reference to the *m_decayData* member. */

        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
        virtual std::string toXMLListExtraAttributes( void ) const ;
        virtual void toXMLListExtraElements( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;

        friend class Unorthodox;
};

/*
============================================================
======================== GaugeBoson ========================
============================================================
*/

class GaugeBoson : public Particle {

    public:
        GaugeBoson( HAPI::Node const &a_node, Database *a_DB, Database *a_parent );
        virtual ~GaugeBoson( );
};

/*
============================================================
========================== Lepton ==========================
============================================================
*/

class Lepton : public Particle {

    private:
        std::string m_generation;                                               /**< The generation of the lepton (i.e., electronic, muonic or tauonic). */

    public:
        Lepton( HAPI::Node const &a_node, Database *a_DB, Database *a_parent );
        virtual ~Lepton( );

        std::string const &generation( void ) const { return( m_generation ); }
        virtual std::string toXMLListExtraAttributes( void ) const ;
};

/*
============================================================
========================== Baryon ==========================
============================================================
*/

class Baryon : public Particle {

    public:
        Baryon( HAPI::Node const &a_node, Database *a_DB, Database *a_parent );
        virtual ~Baryon( );
};

/*
============================================================
======================== Unorthodox ========================
============================================================
*/

class Unorthodox : public Particle {

    public:
        Unorthodox( HAPI::Node const &a_node, Database *a_DB, Database *a_parent );
        virtual ~Unorthodox( );
};

/*
============================================================
========================== Nucleus =========================
============================================================
*/

class Nucleus : public Particle {

    private:
        Nuclide *m_nuclide;                                 /**< The parent nuclide of *this*. */
        int m_Z;                                            /**< The atomic number of the parent nuclide. */
        int m_A;                                            /**< The atomic mass number of the parent nuclide. */
        std::string m_levelName;                            /**< The string representationn of *m_levelIndex*. */
        int m_levelIndex;                                   /**< The index of the excited nucleus state. */
        PQ_suite m_energy;                                  /**< A suite storing the physical quantities representing the nucleus excited energy for the particle.. */

    public:
        Nucleus( HAPI::Node const &node, Database *a_DB, Nuclide *a_parent );
        virtual ~Nucleus( );

        Nuclide const *nuclide( ) const { return( m_nuclide ); }                    /**< Returns a *const* reference to the *m_nuclide* member. */
        int Z( void ) const { return( m_Z ); }                                      /**< Returns a *const* reference to the *m_Z* member. */
        int A( void ) const { return( m_A ); }                                      /**< Returns a *const* reference to the *m_A* member. */
        std::string const &levelName( ) const { return( m_levelName ); }            /**< Returns a *const* reference to the *m_levelName* member of *this*. */
        int levelIndex( void ) const { return( m_levelIndex ); }                    /**< Returns a *const* reference to the *m_levelIndex* member. */
        std::string const &atomsID( void ) const ;

        double massValue( char const *a_unit ) const ;
        PQ_suite const &energy( void ) const { return( m_energy ); }                /**< Returns a *const* reference to the *m_energy* member. */
        double energy( std::string const &a_unit ) const ;
        virtual std::string toXMLListExtraAttributes( void ) const ;
        virtual void toXMLListExtraElements( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
========================== Nuclide =========================
============================================================
*/

class Nuclide : public Particle {

    private:
        Isotope *m_isotope;                                 /**< A pointer to the parent isotope. */
        Nucleus m_nucleus;                                  /**< The nucleus for *this* nuclide. */
        GammaDecayData m_gammaDecayData;                    /**< . */

    public:
        Nuclide( HAPI::Node const &a_node, Database *a_DB, Isotope *a_parent );
        virtual ~Nuclide( );

        int Z( void ) const;
        int A( void ) const;
        std::string const &levelName( void ) const { return( m_nucleus.levelName( ) ); }
                                                                                /**< Returns the result of calling m_nucleus.levelName( ). */
        int levelIndex( void ) const { return( m_nucleus.levelIndex( ) ); }     /**< Returns the result of calling m_nucleus.levelIndex( ). */
        std::string const &atomsID( ) const ;
        std::string const &kind( ) const { return( m_gammaDecayData.kind( ) ); }
        GammaDecayData const &gammaDecayData( ) const { return( m_gammaDecayData ); }

        Isotope const *isotope( ) const { return( m_isotope ); }                /**< Returns a *const* reference to the *m_isotope* member. */
        Nucleus const &nucleus( ) const { return( m_nucleus ); }                /**< Returns a *const* reference to the *m_nucleus* member. */

        PQ_suite const &baseMass( void ) const ;
        double massValue( char const *a_unit ) const ;
        double levelEnergy( std::string const &a_unit ) const { return( m_nucleus.energy( a_unit ) ); }
                                                                                /**< Returns the result of calling m_nucleus.energy( a_unit ). */

        void calculateNuclideGammaBranchStateInfos( PoPI::Database const &a_pops, NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos,
                bool a_alwaysAdd = false ) const ;
        virtual void toXMLListExtraElements( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
========================= Isotope ==========================
============================================================
*/

class Isotope : public SymbolBase {

    private:
        ChemicalElement *m_chemicalElement;         /**< A pointer to the parent chemical element. */
        int m_Z;                                    /**< A atomic number for the isotope. */
        int m_A;                                    /**< The atomic mass number for the isotope. */
        Suite<Nuclide, Isotope> m_nuclides;         /**< The suite of nuclides for this isotope. */

    public:
        Isotope( HAPI::Node const &a_node, Database *a_DB, ChemicalElement *a_parent );
        virtual ~Isotope( );

        ChemicalElement const *chemicalElement( ) const { return( m_chemicalElement ); }    /**< Returns a *const* reference to the *m_isotope* member. */
        int Z( void ) const { return( m_Z ); }                                  /**< Returns the value of the *m_Z* member. */
        int A( void ) const { return( m_A ); }                                  /**< Returns the value of the *m_A* member. */
        Suite<Nuclide, Isotope> const &nuclides( ) const { return( m_nuclides ); }          /**< Returns a *const* reference to the *m_nuclides* member. */

        void calculateNuclideGammaBranchStateInfos( PoPI::Database const &a_pops, NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos ) const ;
        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
===================== ChemicalElement ======================
============================================================
*/

class ChemicalElement : public SymbolBase {

    private:
        int m_Z;                                        /**< A atomic number for all isotopes in *thie* chemical element. */
        std::string m_name;                             /**< The name of the chemical element. */
        Suite<Isotope, ChemicalElement> m_isotopes;     /**< The suite of isotopes for this chemical element. */

    public:
       ChemicalElement( HAPI::Node const &a_node, Database *a_DB, Database *a_parent );
        virtual ~ChemicalElement( );

        int Z( void ) const { return( m_Z ); }                          /**< Returns the value of the *m_Z* member. */
        std::string const &name( void ) const { return( m_name ); }     /**< Returns the value of the *m_name* member. */

        Suite<Isotope, ChemicalElement> const &isotopes( ) const { return( m_isotopes ); }  /**< Returns a *const* reference to the *m_isotopes* member. */

        void calculateNuclideGammaBranchStateInfos( PoPI::Database const &a_pops, NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos ) const ;
        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
=========================== Alias ==========================        // FIXME, there should be an alias base class that Alias and MetaStable inherit from.
============================================================
*/

class Alias : public IDBase {

    private:
        std::string m_pid;                      /**< The id of the particle *this* is an alias for. */
        std::size_t m_pidIndex;                         /**< The index of the particle with id *m_pid*. */

    public:
        Alias( HAPI::Node const &a_node, Database *a_DB, Particle_class a_class = Particle_class::alias );
        virtual ~Alias( );

        std::string const &pid( void ) const { return( m_pid ); }       /**< Returns a *const* reference to the *m_pid* member of *this*. */
        std::size_t pidIndex( void ) const { return( m_pidIndex ); }            /**< Returns a *const* reference to the *m_pidIndex* member of *this*. */
        void setPidIndex( std::size_t a_index ) { m_pidIndex = a_index; }       /**< Set the member *m_pidIndex* to *a_index*. */

        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
======================== MetaStable ========================
============================================================
*/

class MetaStable : public Alias {

    private:
        int m_metaStableIndex;                  /**< The meta-stable index for *this*. */

    public:
        MetaStable( HAPI::Node const &a_node, Database *a_DB );
        virtual ~MetaStable( );

        int metaStableIndex( void ) const { return( m_metaStableIndex ); }                          /**< Returns the value of the *m_metaStableIndex* member. */
        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
};

/*
============================================================
========================= Database =========================
============================================================
*/

class Database {

    private:
        LUPI::FormatVersion m_formatVersion;                            /**< The **GNDS** **format** attribute of the first file read. */
        std::string m_name;                                             /**< The **GNDS** **name** of the first file read in. */
        std::string m_version;                                          /**< The **GNDS** **version** of the first file read in. */
        ParticleList m_list;                                            /**< The internal list of the particles. */
        std::map<std::string, std::size_t> m_idsMap;    // Be careful with this as a map[key] will add key if it is not in the map.
        std::map<int, std::size_t> m_intidsMap;         // Be careful with this as a map[key] will add key if it is not in the map.
                                                                        /**< This maps each particle id to a unique index. */
        SymbolList m_symbolList;                                        /**< The internal list of the symbols. */
        std::map<std::string, std::size_t> m_symbolMap; // Be careful with this as a map[key] will add key if it is not in the map.
                                                                        /**< This maps each symbol to a unique index. */

        std::vector<Alias *> m_unresolvedAliases;                       /**< This is used internally to store aliases when a **PoPs** node is being parsed as the aliases onde is parsed before the particles are parsed. */
        std::vector<Alias *> m_aliases;                                 /**< Represents the **PoPs** *aliases* node which contains a list of the **PoPs** **alias** and **metaStable** nodes. */

        Suite<GaugeBoson, Database> m_gaugeBosons;                      /**< Represents the **PoPs** **gaugeBosons** node which contains a list of the **PoPs** **gaugeBoson** nodes. */
        Suite<Lepton, Database> m_leptons;                              /**< Represents the **PoPs** **leptons** node which contains a list of the **PoPs** **lepton** nodes. */
        Suite<Baryon, Database> m_baryons;                              /**< Represents the **PoPs** **baryons** node which contains a list of the **PoPs** **baryon** nodes. */
        Suite<ChemicalElement, Database> m_chemicalElements;            /**< Represents the **PoPs** **chemicalElements** node which contains a list of the **PoPs** **chemicalElement** nodes. */
        Suite<Unorthodox, Database> m_unorthodoxes;                     /**< Represents the **PoPs** **unorthodoxes** node which contains a list of the **PoPs** **unorthodox** nodes. */

    public:
        Database( );
        Database( std::string const &a_fileName );
        Database( HAPI::Node const &a_database );
        ~Database( );

        LUPI::FormatVersion const &formatVersion( void ) const { return( m_formatVersion ); }   /**< Returns a *const* *reference* to the *m_formatVersion* variable of *this*. */
        std::string const &name( void ) const { return( m_name ); }                         /**< Returns a *const* *reference* to the *m_name* variable of *this*. */
        std::string const &version( void ) const { return( m_version ); }                   /**< Returns a *const* *reference* to the *m_version* variable of *this*. */

        std::vector<Alias *> unresolvedAliases( ) { return( m_unresolvedAliases ); }      /**< Returns a *reference* to the *m_unresolvedAliases* member of *this*. */
        std::size_t numberOfUnresolvedAliases( ) { return( m_unresolvedAliases.size( ) ); }     /**< Returns the number of unresolved aliases. */
        std::vector<std::string> unresolvedAliasIds( ) const ;
        std::vector<Alias *> &aliases( ) { return( m_aliases ); }                            /**< Returns a *const* *reference* to the *m_aliases* variable of *this*. */

        void addFile( char const *a_fileName, bool a_warnIfDuplicate );
        void addFile( std::string const &a_fileName, bool a_warnIfDuplicate );
        void addDatabase( std::string const &a_string, bool a_warnIfDuplicate );
        void addDatabase( HAPI::Node const &a_database, bool a_warnIfDuplicate );
        void addAlias( Alias *a_alias ) { m_aliases.push_back( a_alias ); }                 /**< Added the **Alias** *a_alias* to *this*. */

        std::string::size_type size( void ) const { return( m_list.size( ) ); }             /**< Returns the number of particle in *this*. */
        ParticleList const &list( ) { return( m_list ); }                                    /**< Returns a *const* *reference* to the *m_list* member. */
        SymbolList const &symbolList( ) { return( m_symbolList ); }                          /**< Returns a *const* *reference* to the *m_symbolList* member. */
        std::size_t operator[]( std::string const &a_id ) const ;
        template<typename T> T const &get( std::string const &a_id ) const ;
        template<typename T> T const &get( std::size_t a_index ) const ;
        Particle const &particle( std::string const &a_id ) const { return( get<Particle>( a_id ) ); }  /**< Returns a *const* *reference* to the particle with id *a_id*. */
        Particle const &particle( std::size_t a_index ) const { return( get<Particle>( a_index ) ); }   /**< Returns a *const* *reference* to the particle with index *a_index*. */
        IDBase const &idBase( std::string const &a_id ) const { return( get<IDBase>( a_id ) ); }        /**< Returns a *const* *reference* to a **IDBase** instance with id *a_id*. */
        IDBase const &idBase( std::size_t &a_index ) const { return( get<IDBase>( a_index ) ); }        /**< Returns a *const* *reference* to a **IDBase** instance with id *a_index*. */
        ParticleList const &particleList( ) const { return( m_list ); }                                 /**< Returns a *const* *reference* to the *m_list* variable of *this*. */
        SymbolList symbolList( ) const { return( m_symbolList ); }                                      /**< Returns a *const* *reference* to the *m_symbolList* variable of *this*. */

        bool exists( std::string const &a_id ) const ;
        bool exists( std::size_t a_index ) const ;
        bool existsIntid( int a_intid ) const ;

        Suite<ChemicalElement, Database> const &chemicalElements( ) const { return( m_chemicalElements ); }
                                                                                            /**< Returns a *const* *reference* to the *m_chemicalElements* variable of *this*. */

        bool isParticle( std::string const &a_id ) const { return( get<Base>( a_id ).isParticle( ) ); } /**< Returns **true** if *a_id* is a particle and **false** otherwise. */
        bool isParticle( std::size_t a_index ) const { return( m_list[a_index]->isParticle( ) ); }      /**< Returns **true** if *a_index* is a particle and **false** otherwise. */
        bool isAlias( std::string const &a_id ) const { return( get<Base>( a_id ).isAlias( ) ); }       /**< Returns **true** if *a_id* is an alias and **false** otherwise. */
        bool isAlias( std::size_t a_index ) const { return( m_list[a_index]->isAlias( ) ); }            /**< Returns **true** if *a_index* is an alias and **false** otherwise. */
        bool isMetaStableAlias( std::string const &a_id ) const { return( get<Base>( a_id ).isMetaStableAlias( ) ); }
                                                                                                        /**< Returns **true** if *a_id* is a meta-stable and **false** otherwise. */
        bool isMetaStableAlias( std::size_t a_index ) const { return( m_list[a_index]->isMetaStableAlias( ) ); }
                                                                                                        /**< Returns **true** if *a_index* is a meta-stable and **false** otherwise. */
        std::vector<std::string> aliasReferences( std::string const &a_id );

        std::string final( std::string const &a_id, bool a_returnAtMetaStableAlias = false ) const ;
        std::size_t final( std::size_t a_index, bool a_returnAtMetaStableAlias = false ) const ;

        std::string chemicalElementSymbol( std::string const &a_id ) const ;
        std::string isotopeSymbol( std::string const &a_id ) const ;
        int intid( std::string const &a_id ) const ;
        int intid( std::size_t a_index ) const ;
        std::size_t indexFromIntid( int a_intid ) const ;

        std::size_t add( Base *a_item );
        std::size_t addSymbol( SymbolBase *a_item );

        void calculateNuclideGammaBranchStateInfos( NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos, Database const *a_pops2,
                std::vector<std::string> &a_extraGammaBranchStates ) const ;
        void calculateNuclideGammaBranchStateInfos2( NuclideGammaBranchStateInfos &a_nuclideGammaBranchStateInfos ) const ;

        double massValue( std::string const &a_id, std::string const &a_unit ) const ;

        void saveAs( std::string const &a_fileName ) const ;
        void toXMLList( std::vector<std::string> &a_XMLList, std::string const &a_indent1 ) const ;
        void print( bool a_printIndices );
};

/* *********************************************************************************************************//**
 * Returns the partile in *this* that has index *a_index*.
 *
 * @param a_index                       [in]    The index of the particle to return.
 *
 * @return                                      A *const* reference to the particle at index *a_index*.
 ***********************************************************************************************************/

template<typename T> T const &Database::get( std::size_t a_index ) const {

    Base *particle = m_list[a_index];
    if( particle == nullptr ) throw std::range_error( std::string( "particle not in database" ) );
    T const *object = dynamic_cast<T const *>( particle );
    if( object == nullptr ) throw std::bad_cast( );

    return( *object );
}

/* *********************************************************************************************************//**
 * Returns the partile in *this* that has index *a_index*.
 *
 * @param a_id                          [in]    The **PoPs** id of the particle to return.
 *
 * @return                                      A *const* reference to the particle with id *a_id*.
 ***********************************************************************************************************/

template<typename T> T const &Database::get( std::string const &a_id ) const {

    auto index = (*this)[a_id];
    Base *particle = m_list[index];
    T const *object = dynamic_cast<T const *>( particle );
    if( object == nullptr ) throw std::bad_cast( );

    return( *object );
}

double getPhysicalQuantityAsDouble( PhysicalQuantity const &a_physicalQuantity );
double getPhysicalQuantityOfSuiteAsDouble( PQ_suite const &a_suite, bool a_allowEmpty = false, double a_emptyValue = 0.0 );
bool supportedFormat( LUPI::FormatVersion const &a_formatVersion );
std::string baseAntiQualifierFromID( std::string const &a_id, std::string &a_anti, std::string *a_qualifier = nullptr );

int maximumChemicalElementZ( );
std::string chemicalElementInfoFromZ( int a_Z, bool a_wantSymbol, bool a_asNucleus = false );
std::string const &chemicalElementSymbolFromZ( int a_Z );
int Z_FromChemicalElementSymbol( std::string const &a_symbol );

int family2Integer( Particle_class a_family );
int intidHelper( bool a_isAnti, Particle_class a_family, int a_SSSSSSS );

}

#endif      // End of PoPI_hpp_included
