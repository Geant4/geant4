/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef MCGIDI_hpp_included
#define MCGIDI_hpp_included 1

#define MCGIDI_USE_DOUBLES 1
#ifndef MCGIDI_USE_DOUBLES
    #define MCGIDI_FLOAT float
    #define DATA_MEMBER_VECTOR_FLOAT_OR_DOUBLE DATA_MEMBER_VECTOR_FLOAT
#else
    #define MCGIDI_FLOAT double
    #define DATA_MEMBER_VECTOR_FLOAT_OR_DOUBLE DATA_MEMBER_VECTOR_DOUBLE
#endif

#define _USE_MATH_DEFINES
#include "math.h"

#include <LUPI.hpp>
#include <PoPI.hpp>
#include <GIDI.hpp>

#ifdef MCGIDI_USE_VIRTUAL_FUNCTIONS
    #define MCGIDI_VIRTUAL_FUNCTION virtual
    #define MCGIDI_TRUE_VIRTUAL = 0
#else
    #define MCGIDI_VIRTUAL_FUNCTION
    #define MCGIDI_TRUE_VIRTUAL
#endif

#if defined(HAVE_HIP) && defined(__HIP_DEVICE_COMPILE__)
    #define MCGIDI_PRINTF(...)
#else
    #define MCGIDI_PRINTF printf
#endif

namespace MCGIDI {

class Protare;
class ProtareSingle;
class ProtareComposite;
class ProtareTNSL;
class Reaction;
class OutputChannel;
class ACE_URR_probabilityTables;

}           // End of namespace MCGIDI.

#include <LUPI_dataBuffer.hpp>
#include "MCGIDI_sampling.hpp"
#include "MCGIDI_vector.hpp"
#include "MCGIDI_string.hpp"

namespace MCGIDI {

#define MCGIDI_nullReaction -10001

// FIXME, this should not be used once physicalQuantity can handle changing units.
#define MCGIDI_speedOfLight_cm_sh 299.792458
#define MCGIDI_speedOfLight_cm_sec ( MCGIDI_speedOfLight_cm_sh * 1e8 )
#define MCGIDI_classicalElectronRadius 0.2817940322010228 // Classical electron radius in unit of sqrt( b ).

#define MCGIDI_particleBeta( a_mass_unitOfEnergy, a_kineticEnergy ) \
    ( (a_mass_unitOfEnergy) == 0.0 ? 1 : sqrt( (a_kineticEnergy) * ( (a_kineticEnergy) + 2.0 * (a_mass_unitOfEnergy) ) ) / ( (a_kineticEnergy) + (a_mass_unitOfEnergy) ) )

LUPI_HOST_DEVICE double particleKineticEnergy( double a_mass_unitOfEnergy, double a_particleBeta );
LUPI_HOST_DEVICE double particleKineticEnergyFromBeta2( double a_mass_unitOfEnergy, double a_particleBeta2 );    // a_particleBeta2 = a_particleBeta^2.
LUPI_HOST_DEVICE double boostSpeed( double a_massProjectile, double a_kineticEnergyProjectile, double a_massTarget );
LUPI_HOST_DEVICE int muCOM_From_muLab( double a_muLab, double a_boostBeta, double a_productBeta, double &a_muPlus, double &a_JacobianPlus, double &a_muMinus, double &a_JacobianMinus );

enum class ProtareType { single, composite, TNSL };

namespace Transporting {

enum class URR_mode { none, pdfs, ACE_URR_probabilityTables };

namespace LookupMode {

    enum class Data1d { continuousEnergy, multiGroup };
    enum class Distribution { pdf_cdf, epbs };

}           // End of namespace LookupMode.

namespace Reaction {

    enum class Type { Reactions, OrphanProducts };

}           // End of namespace Reaction.

/*
============================================================
============================ MC ============================
============================================================
*/
class MC : public GIDI::Transporting::Settings {

    private:
        GIDI::Styles::Suite const *m_styles;                        /**< FIXME. */
        std::string m_label;                                        /**< FIXME. */
        double m_energyDomainMax;                                   /**< All reactions with a threshold greater than or equal to this value are ignored. */
        bool m_ignoreENDF_MT5;                                      /**< If true, the ENDF MT 5 reaction is ignored. */
        bool m_sampleNonTransportingParticles;                      /**< If true, all products are sampled, otherwise only transporting particles are sampled. */
        bool m_useSlowerContinuousEnergyConversion;                 /**< If true, the old slower conversion of GIDI::ProtareSingle to MCGIDI::HeatedCrossSectionContinuousEnergy is used. */
        bool m_addExpectedValueData;                                /**< If true, exected value data are included in the construction of a HeatedCrossSectionContinuousEnergy class. */
        LookupMode::Data1d m_crossSectionLookupMode;                /**< Determines how cross sections are evaluated. */
        LookupMode::Data1d m_other1dDataLookupMode;                 /**< Determines how 1d data other than cross sections are evaluated. */
        LookupMode::Distribution m_distributionLookupMode;          /**< Determines how distributions are evaluated and sampled. Currently, only pdf_cdf is allowed. */
        Sampling::Upscatter::Model m_upscatterModel;                /**< FIXME. */
        std::string m_upscatterModelALabel;                         /**< FIXME. */
        URR_mode m_URR_mode;                                        /**< Selects if URR data are to be used, and it so, which type. */
        bool m_wantTerrellPromptNeutronDistribution;                /**< If true, prompt fission neutron distributions are sampled from the Terrell mode. */
        bool m_wantRawTNSL_distributionSampling;                    /**< If true, the TNSL neutron distributions for coherent and incoherent elastic scattering are sampled from the double differential data. Otherwise, they are sampled from the distribution data. */
        std::vector<double> m_fixedGridPoints;                      /**< FIXME. */
        bool m_makePhotonEmissionProbabilitiesOne;                  /**< If true, all photon emission probabilities are set to 1.0 (i.e., all ICCs are set to 0.0). Default is false. */
        bool m_zeroNuclearLevelEnergyWidth;                         /**< If true, all NuclideGammaBranchStateInfo.m_nuclearLevelEnergyWidth members are set to 0.0. Default is false. */

    public:
        LUPI_HOST MC( PoPI::Database const &a_pops, std::string const &a_projectileID, GIDI::Styles::Suite const *a_styles, std::string const &a_label, 
                GIDI::Transporting::DelayedNeutrons a_delayedNeutrons, double energyDomainMax );
        LUPI_HOST MC( PoPI::Database const &a_pops, GIDI::Protare const &a_protare, std::string const &a_label,
                GIDI::Transporting::DelayedNeutrons a_delayedNeutrons, double energyDomainMax );

        LUPI_HOST GIDI::Styles::Suite const *styles( ) const { return( m_styles ); }               /**< Returns the value of the **m_styles**. */
        LUPI_HOST void styles( GIDI::Styles::Suite const *a_styles ) { m_styles = a_styles; }      /**< This is needed for ProtareTNSL, but should be avoided otherwise. FIXME, need to have a better way. */

        LUPI_HOST std::string label( ) const { return( m_label ); }

// FIXME (1) should this not be something like
// GIDI::Styles::Suite const &suite( ) const { return( *m_styles ); }                   /**< Returns a reference to **m_styles**. */
        LUPI_HOST double energyDomainMax( ) const { return( m_energyDomainMax ); }                /**< Returns the value of the **m_energyDomainMax**. */

        LUPI_HOST bool ignoreENDF_MT5( ) const { return( m_ignoreENDF_MT5 ); }                    /**< Returns the value of the **m_ignoreENDF_MT5**. */
        LUPI_HOST void set_ignoreENDF_MT5( bool a_ignoreENDF_MT5 ) { m_ignoreENDF_MT5 = a_ignoreENDF_MT5; }
        LUPI_HOST void setIgnoreENDF_MT5( bool a_ignoreENDF_MT5 ) { set_ignoreENDF_MT5( a_ignoreENDF_MT5 ); }
                            /**< This function is deprecated. Use **set_ignoreENDF_MT5** instead. */

        LUPI_HOST bool sampleNonTransportingParticles( ) const { return( m_sampleNonTransportingParticles); }
        LUPI_HOST void sampleNonTransportingParticles( bool a_sampleNonTransportingParticles ) {
                LUPI::deprecatedFunction( "MCGIDI::Transporting::MC::sampleNonTransportingParticles", "MCGIDI::Transporting::MC::setSampleNonTransportingParticles", "" );
                setSampleNonTransportingParticles( a_sampleNonTransportingParticles ); }
        LUPI_HOST void setSampleNonTransportingParticles( bool a_sampleNonTransportingParticles ) { m_sampleNonTransportingParticles = a_sampleNonTransportingParticles; }

        LUPI_HOST bool useSlowerContinuousEnergyConversion( ) const { return( m_useSlowerContinuousEnergyConversion ); }
        LUPI_HOST void setUseSlowerContinuousEnergyConversion( bool a_useSlowerContinuousEnergyConversion ) {
                m_useSlowerContinuousEnergyConversion = a_useSlowerContinuousEnergyConversion; }

        LUPI_HOST bool addExpectedValueData( ) const { return( m_addExpectedValueData ); }      /**< Returns the value of the *m_addExpectedValueData* member. */
        LUPI_HOST void setAddExpectedValueData( bool a_addExpectedValueData ) { m_addExpectedValueData = a_addExpectedValueData; }  /**< Set the *m_addExpectedValueData* member to *a_addExpectedValueData*. */

        LUPI_HOST LookupMode::Data1d crossSectionLookupMode( ) const { return( m_crossSectionLookupMode ); }      /**< Returns the value of the **m_crossSectionLookupMode**. */
        LUPI_HOST void setCrossSectionLookupMode( LookupMode::Data1d a_crossSectionLookupMode );
        LUPI_HOST void crossSectionLookupMode( LookupMode::Data1d a_crossSectionLookupMode ) {
                LUPI::deprecatedFunction( "MCGIDI::Transporting::MC::crossSectionLookupMode", "MCGIDI::Transporting::MC::setCrossSectionLookupMode", "" );
                setCrossSectionLookupMode( a_crossSectionLookupMode ); }                                            /**< See method **setCrossSectionLookupMode**. This method is deprecated. */

        LUPI_HOST LookupMode::Data1d other1dDataLookupMode( ) const { return( m_other1dDataLookupMode ); }        /**< Returns the value of the **m_other1dDataLookupMode**. */
        LUPI_HOST void setOther1dDataLookupMode( LookupMode::Data1d a_other1dDataLookupMode );
        LUPI_HOST void other1dDataLookupMode( LookupMode::Data1d a_other1dDataLookupMode ) {
                LUPI::deprecatedFunction( "MCGIDI::Transporting::MC::other1dDataLookupMode", "MCGIDI::Transporting::MC::setOther1dDataLookupMode", "" );
                setOther1dDataLookupMode( a_other1dDataLookupMode ); }                                              /**< See method **setOther1dDataLookupMode**. This method is deprecated. */

        LUPI_HOST LookupMode::Distribution distributionLookupMode( ) const { return( m_distributionLookupMode ); }  /**< Returns the value of the **m_distributionLookupMode**. */
        LUPI_HOST void setDistributionLookupMode( LookupMode::Distribution a_distributionLookupMode );
        LUPI_HOST void distributionLookupMode( LookupMode::Distribution a_distributionLookupMode ) { 
                LUPI::deprecatedFunction( "MCGIDI::Transporting::MC::distributionLookupMode", "MCGIDI::Transporting::MC::setDistributionLookupMode", "" );
                setDistributionLookupMode( a_distributionLookupMode ); }                                            /**< See method **setDistributionLookupMode**. This method is deprecated. */

        LUPI_HOST Sampling::Upscatter::Model upscatterModel( ) const { return( m_upscatterModel ); }                /**< Returns the value of the **m_upscatterModel**. */
        LUPI_HOST void set_upscatterModelA( std::string const &a_upscatterModelALabel );
        LUPI_HOST void setUpscatterModelA( std::string const &a_upscatterModelALabel ) { set_upscatterModelA( a_upscatterModelALabel ); }
                                                                                                                    /**< See method **set_upscatterModelA**. */
        LUPI_HOST std::string upscatterModelALabel( ) const { return( m_upscatterModelALabel ); }                   /**< Returns the value of the **m_upscatterModelALabel**. */
        LUPI_HOST void setUpscatterModelB( ) { m_upscatterModel = Sampling::Upscatter::Model::B; }                  /**< Set member *m_upscatterModel* to Sampling::Upscatter::Model::B. */
        LUPI_HOST void setUpscatterModelBSnLimits( ) { m_upscatterModel = Sampling::Upscatter::Model::BSnLimits; }  /**< Set member *m_upscatterModel* to Sampling::Upscatter::Model::BSnLimits. */
        LUPI_HOST void setUpscatterModelDBRC( ) { m_upscatterModel = Sampling::Upscatter::Model::DBRC; }            /**< Set member *m_upscatterModel* to Sampling::Upscatter::Model::DBRC. */

        LUPI_HOST bool want_URR_probabilityTables( ) const {
                LUPI::deprecatedFunction( "MCGIDI::Transporting::MC::want_URR_probabilityTables", "MCGIDI::Transporting::MC::_URR_mode", "" );
                return( m_URR_mode != URR_mode::none ); }            /**< Returns *false* if *m_URR_mode* is **URR_mode::none** and *true* otherwise. This method is deprecated. Please use **_URR_mode** instead. */
        LUPI_HOST void want_URR_probabilityTables( bool a_want_URR_probabilityTables ) {
                LUPI::deprecatedFunction( "MCGIDI::Transporting::MC::want_URR_probabilityTables", "MCGIDI::Transporting::MC::setURR_mode", "" );
                m_URR_mode = URR_mode::none;
                if( a_want_URR_probabilityTables ) m_URR_mode = URR_mode::pdfs;
        }           /**< If *a_want_URR_probabilityTables* is *true* sets *m_URR_mode* to **URR_mode::pdfs**, otherwise set it to **URR_mode::none**. This method is deprecated. Please use **setURR_mode** instead. */

        LUPI_HOST URR_mode _URR_mode( ) const { return( m_URR_mode ); }                                           /**< Returns the value of the **m_URR_mode** member. */
        LUPI_HOST void setURR_mode( URR_mode a_URR_mode ) { m_URR_mode = a_URR_mode; }                            /**< This methods sets member *m_URR_mode* to *a_URR_mode*. */

        LUPI_HOST bool wantTerrellPromptNeutronDistribution( ) const { return( m_wantTerrellPromptNeutronDistribution ); }
                                    /**< Returns the value of the **m_wantTerrellPromptNeutronDistribution** member. */
        LUPI_HOST void setWantTerrellPromptNeutronDistribution( bool a_wantTerrellPromptNeutronDistribution ) {
                m_wantTerrellPromptNeutronDistribution = a_wantTerrellPromptNeutronDistribution;
        }       /**< Set the *m_wantTerrellPromptNeutronDistribution* member to *a_wantTerrellPromptNeutronDistribution*. */
        LUPI_HOST void wantTerrellPromptNeutronDistribution( bool a_wantTerrellPromptNeutronDistribution ) {
                LUPI::deprecatedFunction( "MCGIDI::Transporting::MC::wantTerrellPromptNeutronDistribution", "MCGIDI::Transporting::MC::setWantTerrellPromptNeutronDistribution", "" );
                setWantTerrellPromptNeutronDistribution( a_wantTerrellPromptNeutronDistribution );
        }       /**< See method *setWantTerrellPromptNeutronDistribution*. This method is deprecated. */

        LUPI_HOST bool wantRawTNSL_distributionSampling( ) const { return( m_wantRawTNSL_distributionSampling ); }
        LUPI_HOST bool wantRawTNSL_distributionSampling( ) { return( m_wantRawTNSL_distributionSampling ); }
        LUPI_HOST void set_wantRawTNSL_distributionSampling( bool a_wantRawTNSL_distributionSampling ) { 
                m_wantRawTNSL_distributionSampling = a_wantRawTNSL_distributionSampling; }

        LUPI_HOST std::vector<double> fixedGridPoints( ) const { return( m_fixedGridPoints ); }
        LUPI_HOST void fixedGridPoints( std::vector<double> a_fixedGridPoints ) { m_fixedGridPoints = a_fixedGridPoints; }

        LUPI_HOST bool makePhotonEmissionProbabilitiesOne( ) const { return( m_makePhotonEmissionProbabilitiesOne ); }
                                            /**< Returns the value of the **m_makePhotonEmissionProbabilitiesOne** member. */
        LUPI_HOST void setMakePhotonEmissionProbabilitiesOne( bool a_makePhotonEmissionProbabilitiesOne ) { m_makePhotonEmissionProbabilitiesOne = a_makePhotonEmissionProbabilitiesOne; }
                                            /**< Sets member *m_makePhotonEmissionProbabilitiesOne* to *a_makePhotonEmissionProbabilitiesOne*. */
        LUPI_HOST bool zeroNuclearLevelEnergyWidth( ) const { return( m_zeroNuclearLevelEnergyWidth ); }
                                            /**< Returns the value of the **m_zeroNuclearLevelEnergyWidth** member. */
        LUPI_HOST void setZeroNuclearLevelEnergyWidth( bool a_zeroNuclearLevelEnergyWidth ) { m_zeroNuclearLevelEnergyWidth = a_zeroNuclearLevelEnergyWidth ; }
                                            /**< Sets member *m_zeroNuclearLevelEnergyWidth* to *a_zeroNuclearLevelEnergyWidth*. */

        LUPI_HOST void process( GIDI::Protare const &a_protare );
};

}           // End of namespace Transporting.

enum class TwoBodyOrder { notApplicable, firstParticle, secondParticle };

/*
============================================================
============= ACE_URR_probabilityTablesFromGIDI ============
============================================================
*/

class ACE_URR_probabilityTablesFromGIDI {

    public:
        LUPI_HOST ACE_URR_probabilityTablesFromGIDI( );
        LUPI_HOST ~ACE_URR_probabilityTablesFromGIDI( );

        std::map<std::string, ACE_URR_probabilityTables *> m_ACE_URR_probabilityTables;     // The string is the reaction's label.
};

/*
============================================================
========================= SetupInfo ========================
============================================================
*/
class SetupInfo {

    public:
        ProtareSingle &m_protare;                               /**< The protare the data are loaded into. */
        GIDI::ProtareSingle const &m_GIDI_protare;              /**< The GIDI protare the data are loaded from. */
        PoPI::Database const &m_popsUser;                       /**< The PoPs from the user. */
        PoPI::Database const &m_pops;                           /**< The PoPs from the GNDS protare. */
        int m_neutronIndex;
        int m_photonIndex;
        LUPI::FormatVersion m_formatVersion;
        double m_Q;
        double m_product1Mass;
        double m_product2Mass;
        double m_domainMin;
        double m_domainMax;
        TwoBodyOrder m_twoBodyOrder;
        bool m_isPairProduction;
        bool m_isPhotoAtomicIncoherentScattering;
        std::string m_distributionLabel;                        /**< Set by the ProtareSingle constructor to the distribution label to use for all products. */
        std::map<std::string, int> m_particleIntids;            /**< A list of the particle intids for the transportable particles. */
        std::map<std::string, int> m_particleIndices;           /**< A list of the particle indices for the transportable particles. */
        GIDI::Reaction const *m_reaction;                       /**< A pointer to the current reaction whose data are being filled from the GIDI::ProtareSingle reaction. */
        Transporting::Reaction::Type m_reactionType;
        int m_initialStateIndex;                                /**< If not -1, then reaction contains a branching gamma data with this index for the data in m_nuclideGammaBranchStateInfos member of its **ProtareSingle** instance. */
        bool m_hasFinalStatePhotons;                            /**< If **true**, the reaction has a photon with finalState attribute. */
        std::map<std::string, int> m_initialStateIndices;       /**< If not -1, then reaction contains a branching gamma data with this index for the data in m_nuclideGammaBranchStateInfos member of its **ProtareSingle** instance. */
        std::map<std::string, int> m_stateNamesToIndices;       /**< A map of nuclide PoPs ids to their index in the ProtareSingle::m_nuclideGammaBranchStateInfos member. */
        std::map<std::string, double> m_nuclearLevelEnergies;   /**< A map of nuclide PoPs id to their nuclear level energy. */
        std::map<std::string, ACE_URR_probabilityTablesFromGIDI *> m_ACE_URR_probabilityTablesFromGIDI;
        GIDI::GRIN::GRIN_continuumGammas const *m_GRIN_continuumGammas;

        LUPI_HOST SetupInfo( ProtareSingle &a_protare, GIDI::ProtareSingle const &a_GIDI_protare, PoPI::Database const &a_popsUser, 
                PoPI::Database const &a_pops );
        LUPI_HOST ~SetupInfo( );
};

/*
============================================================
=========================== Others =========================
============================================================
*/
LUPI_HOST int MCGIDI_popsIntid( PoPI::Database const &a_pops, std::string const &a_ID );
LUPI_HOST int MCGIDI_popsIndex( PoPI::Database const &a_pops, std::string const &a_ID );

#if 0
/* *********************************************************************************************************//**
 * This function does a binary search of *a_Xs* for the *index* for which *a_Xs*[*index*] <= *a_x* < *a_Xs*[*index*+1].
 * The values of *a_Xs* must be ascending (i.e., *a_Xs*[i] < *a_Xs*[i+1]).
 *
 *
 *   Returns -2 if a_x < a_Xs[0] or 0 if a_boundIndex is true,
 *           -1 if a_x > last point of a_Xs or a_Xs.size( ) - 1 if a_boundIndex is true, or
 *           the lower index of a_Xs which bound a_x otherwise.
 *
 * Note, when *a_boundIndex* is false the returned *index* can be negative and when it is true
 * the return value will be a valid index of *a_Xs*, including its last point. The index of the last
 * point is only returned when *a_boundIndex* is true and *a_x* is great than the last point of *a_Xs*.
 *
 * @param a_x               [in]    The values whose bounding index within *a_Xs* is to be determined.
 * @param a_Xs              [in]    The list of ascending values.
 * @param a_boundIndex      [in]    If true, out-of-bounds values a treated as end points.
 *
 * @return                          The *index*.
 ***********************************************************************************************************/
#endif

LUPI_HOST_DEVICE inline int binarySearchVector( double a_x, Vector<double> const &a_Xs, bool a_boundIndex = false ) {

    int lower = 0, middle, upper = (int) a_Xs.size( ) - 1;

    if( a_x < a_Xs[0] ) {
        if( a_boundIndex ) return( 0 );
        return( -2 );
    }

    if( a_x > a_Xs[upper] ) {
        if( a_boundIndex ) return( upper );
        return( -1 );
    }

    while( 1 ) {
        middle = ( lower + upper ) >> 1;
        if( middle == lower ) break;
        if( a_x < a_Xs[middle] ) {
            upper = middle; }
        else {
            lower = middle;
        }
    }
    return( lower );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE inline int binarySearchVectorBounded( double a_x, Vector<double> const &a_Xs, int a_lower, 
                int a_upper, bool a_boundIndex ) {

    int middle;

    if( a_x < a_Xs[a_lower] ) {
        if( a_boundIndex ) return( 0 );
        return( -2 );
    }

    if( a_x > a_Xs[a_upper] ) {
        if( a_boundIndex ) return( a_upper );
        return( -1 );
    }

    while( 1 ) {
        middle = ( a_lower + a_upper ) >> 1;
        if( middle == a_lower ) break;
        if( a_x < a_Xs[middle] ) {
            a_upper = middle; }
        else {
            a_lower = middle;
        }
    }
    return( a_lower );
}

}           // End of namespace MCGIDI.

#include "MCGIDI_functions.hpp"
#include "MCGIDI_distributions.hpp"

namespace MCGIDI {

enum class ChannelType { none, twoBody, uncorrelatedBodies };

/*
============================================================
===================== MultiGroupHash =======================
============================================================
*/
class MultiGroupHash {

    private:
        Vector<double> m_boundaries;                                    /**< The list of multi-group boundaries. */

        LUPI_HOST void initialize( GIDI::Protare const &a_protare, GIDI::Styles::TemperatureInfo const &a_temperatureInfo, std::string a_particleID );

    public:
        LUPI_HOST MultiGroupHash( std::vector<double> a_boundaries );
        LUPI_HOST MultiGroupHash( GIDI::Protare const &a_protare, GIDI::Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_particleID = "" );
        LUPI_HOST MultiGroupHash( GIDI::Protare const &a_protare, GIDI::Transporting::Particles const &a_particles );

        LUPI_HOST_DEVICE Vector<double> const &boundaries( ) const { return( m_boundaries ); }   /**< Returns a reference to **m_styles**. */
        LUPI_HOST_DEVICE int index( double a_domain ) const {
            int _index = binarySearchVector( a_domain, m_boundaries );

            if( _index == -2 ) return( 0 );
            if( _index == -1 ) return( m_boundaries.size( ) - 2 );
            return( _index );
        }
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
====================== URR_protareInfo =====================
============================================================
*/
class URR_protareInfo {

    public:
        bool m_inURR;
        double m_rng_Value;

        LUPI_HOST_DEVICE URR_protareInfo( ) : m_inURR( false ), m_rng_Value( 0.0 ) { }
        LUPI_HOST_DEVICE URR_protareInfo( URR_protareInfo const &a_URR_protareInfo ) {
            m_inURR = a_URR_protareInfo.m_inURR;
            m_rng_Value = a_URR_protareInfo.m_rng_Value;
        }
        LUPI_HOST_DEVICE URR_protareInfo &operator=( URR_protareInfo const &a_rhs ) {

            if( this != &a_rhs ) {
                m_inURR = a_rhs.inURR( );
                m_rng_Value = a_rhs.rng_Value( );
            }

            return( *this );
        }

        LUPI_HOST_DEVICE bool inURR( ) const { return( m_inURR ); }
        LUPI_HOST_DEVICE double rng_Value( ) const { return( m_rng_Value ); }
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
===================== URR_protareInfos =====================
============================================================
*/
class URR_protareInfos {

    private:
        Vector<URR_protareInfo> m_URR_protareInfos;

    public:
        LUPI_HOST_DEVICE URR_protareInfos( ) : m_URR_protareInfos( ) { }
        LUPI_HOST URR_protareInfos( Vector<Protare *> &a_protares );

        LUPI_HOST void setup( Vector<Protare *> &a_protares );

        LUPI_HOST_DEVICE std::size_t size( ) const { return( m_URR_protareInfos.size( ) ); }
        LUPI_HOST_DEVICE URR_protareInfo const &operator[]( std::size_t a_index ) const { return( m_URR_protareInfos[a_index] ); }  /**< Returns the instance of *m_URR_protareInfos* at index *a_index*. */
template <typename RNG>
        inline LUPI_HOST_DEVICE void updateProtare( MCGIDI::Protare const *a_protare, double a_energy, RNG && a_rng );

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST_DEVICE long internalSize( ) const { return m_URR_protareInfos.internalSize( ); }
};

/*
============================================================
================= ACE_URR_probabilityTable =================
============================================================
*/

class ACE_URR_probabilityTable {

    public:
        double m_energy;                                                /**< The projectile energy where the data are specified. */
        Vector<double> m_propabilities;                                 /**< The probability for each cross section. */
        Vector<double> m_crossSections;                                 /**< The cross section for each probability. */

        LUPI_HOST_DEVICE ACE_URR_probabilityTable( );
        LUPI_HOST ACE_URR_probabilityTable( double a_energy, std::vector<double> const &a_propabilities, std::vector<double> const &a_crossSection );
        LUPI_HOST_DEVICE ~ACE_URR_probabilityTable( );

        LUPI_HOST_DEVICE double energy( ) const { return( m_energy ); }
        LUPI_HOST_DEVICE double sample( double a_rng_Value );
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
================= ACE_URR_probabilityTables ================
============================================================
*/

class ACE_URR_probabilityTables {

    public:
        Vector<double> m_energies;                                          /**< List of energies where probabilities tables are given. */
        Vector<ACE_URR_probabilityTable *> m_ACE_URR_probabilityTables;     /**< List of probabilities tables. One for each energy in *m_energies*. */

        LUPI_HOST_DEVICE ACE_URR_probabilityTables( );
        LUPI_HOST_DEVICE ACE_URR_probabilityTables( std::size_t a_capacity );
        LUPI_HOST_DEVICE ~ACE_URR_probabilityTables( );

        LUPI_HOST_DEVICE std::size_t capacity( ) const { return( m_energies.capacity( ) ); }
                                                                            /**< Returns the number of energies allocated to store probability tables. */
        LUPI_HOST_DEVICE std::size_t size( ) const { return( m_energies.size( ) ); }
                                                                            /**< Returns the number of energies that have URR probability tables. */
        LUPI_HOST_DEVICE void reserve( std::size_t a_capacity );
        LUPI_HOST_DEVICE void push_back( ACE_URR_probabilityTable *a_ACE_URR_probabilityTable );

        LUPI_HOST_DEVICE double domainMin( ) const { return( m_energies[0] ); }         /**< Returns the minimum energy where URR data are specified. */
        LUPI_HOST_DEVICE double domainMax( ) const { return( m_energies.back( ) ); }    /**< Returns the maximum energy where URR data are specified. */
        LUPI_HOST_DEVICE double sample( double a_energy, double a_rng_Value );

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======== HeatedReactionCrossSectionContinuousEnergy ========
============================================================
*/

class HeatedReactionCrossSectionContinuousEnergy {

    private:
        int m_offset;                                               /**< The offset relative to the cross section grid of the first cross section value in *m_crossSections*. */
        double m_threshold;                                         /**< The threshold for the reaction. */
        Vector<MCGIDI_FLOAT> m_crossSections;                       /**< The reaction's cross section. */
        Transporting::URR_mode m_URR_mode;                          /**< The URR data (i.e., mode) *this* has. */
        Probabilities::ProbabilityBase2d *m_URR_probabilityTables;  /**< Pointer to pdf URR probabilities if they were loaded. */
        ACE_URR_probabilityTables *m_ACE_URR_probabilityTables;     /**< The ACE URR probability tables for the reaction's cross section, if they were loaded. */

    public:
        LUPI_HOST_DEVICE HeatedReactionCrossSectionContinuousEnergy( );
        LUPI_HOST HeatedReactionCrossSectionContinuousEnergy( int a_offset, double a_threshold, Vector<double> &a_crossSection );
        LUPI_HOST HeatedReactionCrossSectionContinuousEnergy( double a_threshold, GIDI::Functions::Ys1d const &a_crossSection, 
                        Probabilities::ProbabilityBase2d *a_URR_probabilityTables, ACE_URR_probabilityTables *a_ACE_URR_probabilityTables );
        LUPI_HOST_DEVICE ~HeatedReactionCrossSectionContinuousEnergy( );

        LUPI_HOST_DEVICE double threshold( ) const { return( m_threshold ); }                           /**< Returns the value of the **m_threshold**. */
        LUPI_HOST_DEVICE int offset( ) const { return( m_offset ); }                                    /**< Returns the value of the **m_offset**. */
        LUPI_HOST Vector<MCGIDI_FLOAT> const &crossSections( ) const { return( m_crossSections ); }     /**< Returns a reference to the member **m_crossSections**. */
        LUPI_HOST_DEVICE bool hasURR_probabilityTables( ) const {
            return( ( m_URR_probabilityTables != nullptr ) || ( m_ACE_URR_probabilityTables != nullptr ) );
        }                                                           /**< Returns true if URR probability tables data present and false otherwise. */
        LUPI_HOST_DEVICE Transporting::URR_mode URR_mode( ) const { return( m_URR_mode ); }             /**< Returns the value of **m_URR_mode**. */
        LUPI_HOST_DEVICE double URR_domainMin( ) const ;
        LUPI_HOST_DEVICE double URR_domainMax( ) const ;
        LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d *URR_probabilityTables( ) const { return( m_URR_probabilityTables ); }      /**< Returns the value of the *m_URR_probabilityTables*. */
        LUPI_HOST_DEVICE ACE_URR_probabilityTables *_ACE_URR_probabilityTables( ) const { return( m_ACE_URR_probabilityTables ); }    /**< Returns the value of the *m_ACE_URR_probabilityTables*. */
        LUPI_HOST_DEVICE double crossSection( std::size_t a_index ) const {
            int index = static_cast<int>( a_index ) - m_offset;
            if( index < 0 ) return( 0.0 );
            if( index >= static_cast<int>( m_crossSections.size( ) ) ) return( 0.0 );

            return( m_crossSections[index] );
        }
        LUPI_HOST GIDI::Functions::XYs1d crossSectionAsGIDI_XYs1d( double a_temperature, Vector<double> const &a_energies ) const ;

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );

        LUPI_HOST void print( ProtareSingle const *a_protareSingle, std::string const &a_indent, std::string const &a_iFormat, 
                std::string const &a_energyFormat, std::string const &a_dFormat ) const ;
};

/*
============================================================
=================== ContinuousEnergyGain ===================
============================================================
*/
class ContinuousEnergyGain {

    private:
        int m_particleIntid;
        int m_particleIndex;
        int m_userParticleIndex;
        Vector<MCGIDI_FLOAT> m_gain;

    public:
        LUPI_HOST_DEVICE ContinuousEnergyGain( );
        LUPI_HOST ContinuousEnergyGain( int a_particleIntid, int a_particleIndex, std::size_t a_size );

        LUPI_HOST ContinuousEnergyGain &operator=( ContinuousEnergyGain const &a_continuousEnergyGain );

        LUPI_HOST_DEVICE int particleIntid( ) const { return( m_particleIntid); }               /**< Returns the value of the *m_particleIntid* member of *this*. */
        LUPI_HOST_DEVICE int particleIndex( ) const { return( m_particleIndex ); }              /**< Returns the value of the *m_particleIndex* member of *this*. */
        LUPI_HOST_DEVICE int userParticleIndex( ) const { return( m_userParticleIndex ); }      /**< Returns the value of the *m_userParticleIndex* member of *this*. */
        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {
                if( a_particleIndex == m_particleIndex ) m_userParticleIndex = a_userParticleIndex; }
                                                        /**< Sets member *m_userParticleIndex* to *a_userParticleIndex* if particle's index matchs *m_particleIndex*. */
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {
                if( a_particleIntid == m_particleIntid ) m_userParticleIndex = a_userParticleIndex; }
                                                        /**< Sets member *m_userParticleIntid* to *a_userParticleIndex* if particle's intid matchs *m_particleIntid*. */
        LUPI_HOST_DEVICE Vector<MCGIDI_FLOAT> const &gain( ) const { return( m_gain ); }
        LUPI_HOST void adjustGain( int a_energy_index, double a_gain ) { m_gain[a_energy_index] += a_gain; }
        LUPI_HOST_DEVICE double gain( int a_energy_index, double a_energy_fraction ) const ;

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST void print( ProtareSingle const *a_protareSingle, std::string const &a_indent, std::string const &a_iFormat, 
                std::string const &a_energyFormat, std::string const &a_dFormat ) const ;
};

/*
============================================================
=========== HeatedCrossSectionContinuousEnergy =============
============================================================
*/
class HeatedCrossSectionContinuousEnergy {

    private:
        double m_temperature;                                   /**< The target temperature of the data. */
        Vector<int> m_hashIndices;                              /**< The indicies for the energy hash function. */
        Vector<double> m_energies;                              /**< Energy grid for cross sections. */
        Vector<MCGIDI_FLOAT> m_totalCrossSection;               /**< The total cross section. */
        Vector<MCGIDI_FLOAT> m_depositionEnergy;                /**< The total continuous energy, deposition-energy cross section (related to the kinetic energy of the untracked outgoing particles). */
        Vector<MCGIDI_FLOAT> m_depositionMomentum;              /**< The total continuous energy, deposition-momentum cross section. */
        Vector<MCGIDI_FLOAT> m_productionEnergy;                /**< The total continuous energy, Q-value cross section. */
        Vector<ContinuousEnergyGain *> m_gains;                 /**< The total continuous energy, gain cross section for each tracked particle. */
        Transporting::URR_mode m_URR_mode;                      /**< The URR data (i.e., mode) *this* has. */
        Vector<int> m_reactionsInURR_region;                    /**< A list of reactions within or below the upper URR regions. This is empty unless URR probability tables present and used. */
        Vector<HeatedReactionCrossSectionContinuousEnergy *> m_reactionCrossSections;
                                                                /**< Reaction cross section data for each reaction. */
        ACE_URR_probabilityTables *m_ACE_URR_probabilityTables; /**< The ACE URR probability tables for the summed URR cross section, if they were loaded. */

    public:
        LUPI_HOST_DEVICE HeatedCrossSectionContinuousEnergy( );
        LUPI_HOST HeatedCrossSectionContinuousEnergy( SetupInfo &a_setupInfo, Transporting::MC const &a_settings, GIDI::Transporting::Particles const &a_particles,
                DomainHash const &a_domainHash, GIDI::Styles::TemperatureInfo const &a_temperatureInfo, std::vector<GIDI::Reaction const *> const &a_reactions,
                std::vector<GIDI::Reaction const *> const &a_orphanProducts, bool a_fixedGrid, bool a_zeroReactions );
        LUPI_HOST_DEVICE ~HeatedCrossSectionContinuousEnergy( );

        LUPI_HOST_DEVICE int evaluationInfo( int a_hashIndex, double a_energy, double *a_energyFraction ) const ;

        LUPI_HOST HeatedReactionCrossSectionContinuousEnergy const *reactionCrossSection( int a_index ) const 
                { return( m_reactionCrossSections[a_index] ); } /**< Returns the reaction cross section at index *a_index*. */

        LUPI_HOST_DEVICE double temperature( ) const { return( m_temperature ); }           /**< Returns the value of the **m_temperature** member. */
        LUPI_HOST_DEVICE double minimumEnergy( ) const { return( m_energies[0] ); }         /**< Returns the minimum cross section domain. */
        LUPI_HOST_DEVICE double maximumEnergy( ) const { return( m_energies.back( ) ); }    /**< Returns the maximum cross section domain. */
        LUPI_HOST_DEVICE int numberOfReactions( ) const { return( (int) m_reactionCrossSections.size( ) ); } 
                                                                /**< Returns the number of reaction cross section. */

        LUPI_HOST_DEVICE int thresholdOffset( int a_reactionIndex ) const { return( m_reactionCrossSections[a_reactionIndex]->offset( ) ); }
                                                                /**< Returns the offset for the cross section for the reaction with index *a_reactionIndex*. */
        LUPI_HOST_DEVICE double threshold( int a_reactionIndex ) const { return( m_reactionCrossSections[a_reactionIndex]->threshold( ) ); }
                                                                /**< Returns the threshold for the reaction with index *a_reactionIndex*. */
        LUPI_HOST_DEVICE bool hasURR_probabilityTables( ) const ;
        LUPI_HOST_DEVICE double URR_domainMin( ) const ;
        LUPI_HOST_DEVICE double URR_domainMax( ) const ;
        LUPI_HOST_DEVICE bool reactionHasURR_probabilityTables( int a_index ) const { return( m_reactionCrossSections[a_index]->hasURR_probabilityTables( ) ); }

        LUPI_HOST_DEVICE Vector<MCGIDI_FLOAT> &totalCrossSection( ) { return( m_totalCrossSection ); }     /**< Returns a reference to member *m_totalCrossSection*. */
        LUPI_HOST_DEVICE double crossSection(                               URR_protareInfos const &a_URR_protareInfos, int a_URR_index, int a_hashIndex, double a_energy, bool a_sampling = false ) const ;
        LUPI_HOST GIDI::Functions::XYs1d crossSectionAsGIDI_XYs1d( ) const ;

        LUPI_HOST_DEVICE double reactionCrossSection(  int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_URR_index, int a_hashIndex, double a_energy, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection2( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_URR_index, double a_energy, int a_energyIndex, double a_energyFraction, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection(  int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_URR_index, double a_energy ) const ;
        LUPI_HOST GIDI::Functions::XYs1d reactionCrossSectionAsGIDI_XYs1d( int a_reactionIndex ) const ;

        LUPI_HOST_DEVICE double depositionEnergy(   int a_hashIndex, double a_energy ) const ;
        LUPI_HOST_DEVICE double depositionMomentum( int a_hashIndex, double a_energy ) const ;
        LUPI_HOST_DEVICE double productionEnergy(   int a_hashIndex, double a_energy ) const ;
        LUPI_HOST_DEVICE double gain(               int a_hashIndex, double a_energy, int a_particleIndex ) const ;
        LUPI_HOST_DEVICE double gainViaIntid(       int a_hashIndex, double a_energy, int a_particleIntid ) const ;

        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );

        LUPI_HOST_DEVICE Vector<double> const &energies( ) const { return( m_energies ); }       /**< Returns a reference to **m_styles**. */

        LUPI_HOST void print( ProtareSingle const *a_protareSingle, std::string const &a_indent, std::string const &a_iFormat, 
                std::string const &a_energyFormat, std::string const &a_dFormat ) const ;
};

/*
============================================================
============ HeatedCrossSectionsContinuousEnergy ===========
============================================================
*/
class HeatedCrossSectionsContinuousEnergy {

    private:
        Vector<double> m_temperatures;                                      /**< The list of temperatures that have **HeatedCrossSectionContinuousEnergy** data. */
        Vector<double> m_thresholds;                                        /**< The threshold for each reaction. */
        Vector<HeatedCrossSectionContinuousEnergy *> m_heatedCrossSections; /**< One **HeatedCrossSectionContinuousEnergy** instance for each temperature in *m_temperature*. */

    public:
        LUPI_HOST_DEVICE HeatedCrossSectionsContinuousEnergy( );
        LUPI_HOST_DEVICE ~HeatedCrossSectionsContinuousEnergy( );

        LUPI_HOST void update( LUPI::StatusMessageReporting &a_smr, SetupInfo &a_setupInfo, Transporting::MC const &a_settings, GIDI::Transporting::Particles const &a_particles, DomainHash const &a_domainHash, 
                GIDI::Styles::TemperatureInfos const &a_temperatureInfos, std::vector<GIDI::Reaction const *> const &a_reactions, 
                std::vector<GIDI::Reaction const *> const &a_orphanProducts, bool a_fixedGrid, bool a_zeroReactions );

        LUPI_HOST_DEVICE double minimumEnergy( ) const { return( m_heatedCrossSections[0]->minimumEnergy( ) ); }
                                                                    /**< Returns the minimum cross section domain. */
        LUPI_HOST_DEVICE double maximumEnergy( ) const { return( m_heatedCrossSections[0]->maximumEnergy( ) ); }
                                                                    /**< Returns the maximum cross section domain. */
        LUPI_HOST_DEVICE Vector<double> const &temperatures( ) const { return( m_temperatures ); }   /**< Returns the value of the **m_temperatures**. */
        Vector<HeatedCrossSectionContinuousEnergy *> &heatedCrossSections( ) { return( m_heatedCrossSections ); }

        LUPI_HOST_DEVICE double threshold( std::size_t a_index ) const { return( m_thresholds[a_index] ); }     /**< Returns the threshold for the reaction at index *a_index*. */
        LUPI_HOST_DEVICE bool hasURR_probabilityTables( ) const { return( m_heatedCrossSections[0]->hasURR_probabilityTables( ) ); }
        LUPI_HOST_DEVICE double URR_domainMin( ) const { return( m_heatedCrossSections[0]->URR_domainMin( ) ); }
        LUPI_HOST_DEVICE double URR_domainMax( ) const { return( m_heatedCrossSections[0]->URR_domainMax( ) ); }
        LUPI_HOST_DEVICE bool reactionHasURR_probabilityTables( int a_index ) const { return( m_heatedCrossSections[0]->reactionHasURR_probabilityTables( a_index ) ); }

        LUPI_HOST_DEVICE double crossSection(                              URR_protareInfos const &a_URR_protareInfos, int a_URR_index, int a_hashIndex, 
                double a_temperature, double a_energy, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE void crossSectionVector( double a_temperature, double a_userFactor, std::size_t a_numberAllocated, 
                double *a_crossSectionVector ) const ;
        LUPI_HOST GIDI::Functions::XYs1d crossSectionAsGIDI_XYs1d( double a_temperature ) const ;

        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_URR_index, int a_hashIndex, 
                double a_temperature, double a_energy, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_URR_index, double a_temperature, double a_energy_in ) const ;
        LUPI_HOST GIDI::Functions::XYs1d reactionCrossSectionAsGIDI_XYs1d( int a_reactionIndex, double a_temperature ) const ;

        template <typename RNG>
        inline LUPI_HOST_DEVICE int sampleReaction(                               URR_protareInfos const &a_URR_protareInfos, int a_URR_index, int a_hashIndex, 
                double a_temperature, double a_energy, double a_crossSection, RNG && a_rng) const ;

        LUPI_HOST_DEVICE double depositionEnergy(   int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double depositionMomentum( int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double productionEnergy(   int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double gain(               int a_hashIndex, double a_temperature, double a_energy, int a_particleIndex ) const ;
        LUPI_HOST_DEVICE double gainViaIntid(       int a_hashIndex, double a_temperature, double a_energy, int a_particleIntid ) const ;

        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );

        LUPI_HOST void print( ProtareSingle const *a_protareSingle, std::string const &a_indent, std::string const &a_iFormat, 
                std::string const &a_energyFormat, std::string const &a_dFormat ) const ;
};

/*
============================================================
====================== MultiGroupGain ======================
============================================================
*/
class MultiGroupGain {

    private:
        int m_particleIntid;
        int m_particleIndex;
        int m_userParticleIndex;
        Vector<double> m_gain;

    public:
        LUPI_HOST_DEVICE MultiGroupGain( );
        LUPI_HOST MultiGroupGain( int a_particleIntid, int a_particleIndex, GIDI::Vector const &a_gain );

        LUPI_HOST MultiGroupGain &operator=( MultiGroupGain const &a_multiGroupGain );

        LUPI_HOST_DEVICE int particleIntid( ) const { return( m_particleIntid ); }              /**< Returns the value of the *m_particleIntid* member of *this*. */
        LUPI_HOST_DEVICE int particleIndex( ) const { return( m_particleIndex ); }              /**< Returns the value of the *m_particleIndex* member of *this*. */
        LUPI_HOST_DEVICE int userParticleIndex( ) const { return( m_userParticleIndex ); }      /**< Returns the value of the *m_userParticleIndex* member of *this*. */
        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex ) {
                if( a_particleIndex == m_particleIndex ) m_userParticleIndex = a_userParticleIndex; }
                                                        /**< Sets member *m_userParticleIndex* to *a_userParticleIndex* if particle's index matchs *m_particleIndex*. */
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex ) {
                if( a_particleIntid == m_particleIntid ) m_userParticleIndex = a_userParticleIndex; }
                                                        /**< Sets member *m_userParticleIntid* to *a_userParticleIndex* if particle's intid matchs *m_particleIntid*. */
        LUPI_HOST_DEVICE Vector<double> const &gain( ) const { return( m_gain ); }
        LUPI_HOST_DEVICE double gain( int a_hashIndex ) const { return( m_gain[a_hashIndex] ); }

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST void write( FILE *a_file ) const ;
};

/*
============================================================
=========== HeatedReactionCrossSectionMultiGroup ===========
============================================================
*/
class HeatedReactionCrossSectionMultiGroup {

    private:
        double m_threshold;
        int m_offset;
        Vector<double> m_crossSections;             // Multi-group reaction cross section
        double m_augmentedThresholdCrossSection;    // Augmented cross section at m_offset for rejecting when projectile energy is below m_threshold.
                                                    // This value is added to m_crossSections[m_offset] when sampling an isotope or reaction.

    public:
        LUPI_HOST_DEVICE HeatedReactionCrossSectionMultiGroup( );
        LUPI_HOST HeatedReactionCrossSectionMultiGroup( SetupInfo &a_setupInfo, Transporting::MC const &a_settings, int a_offset, 
                std::vector<double> const &a_crossSection, double a_threshold );

        LUPI_HOST_DEVICE double operator[]( std::size_t a_index ) const { return( m_crossSections[a_index] ); }  /**< Returns the value of the cross section at multi-group index *a_index*. */
        LUPI_HOST_DEVICE double threshold( ) const { return( m_threshold ); }        /**< Returns the value of the **m_threshold**. */
        LUPI_HOST_DEVICE int offset( ) const { return( m_offset ); }                 /**< Returns the value of the **m_offset**. */
        LUPI_HOST_DEVICE double crossSection( std::size_t a_index, bool a_sampling = false ) const {
            int index = (int)a_index - m_offset;
            if( index < 0 ) return( 0 );
            if( index >= (int)m_crossSections.size( ) ) return( 0 );

            double _crossSection( m_crossSections[index] );
            if( a_sampling && ( index == 0 ) ) {
                _crossSection += m_augmentedThresholdCrossSection;
            }
            return( _crossSection );
        }
        LUPI_HOST_DEVICE double augmentedThresholdCrossSection( ) const { return( m_augmentedThresholdCrossSection ); }  /**< Returns the value of the **m_augmentedThresholdCrossSection**. */
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST void write( FILE *a_file, int a_reactionIndex ) const ;
};

/*
============================================================
============== HeatedCrossSectionMultiGroup  ==============
============================================================
*/
class HeatedCrossSectionMultiGroup {

    private:
        Vector<double> m_totalCrossSection;                 /**< The total multi-group cross section. */
        Vector<double> m_augmentedCrossSection;             /**< The total multi-group cross section used for sampling with rejection (i.e., null-reactions). */
        Vector<double> m_depositionEnergy;                  /**< The total multi-group, deposition-energy cross section (related to the kinetic energy of the untracked outgoing particles). */
        Vector<double> m_depositionMomentum;                /**< The total multi-group, deposition-momentum cross section. */
        Vector<double> m_productionEnergy;                  /**< The total multi-group, Q-value cross section. */
        Vector<MultiGroupGain *> m_gains;                   /**< The total multi-group, gain cross section for each tracked particle. */
        Vector<HeatedReactionCrossSectionMultiGroup *> m_reactionCrossSections;

    public:
        LUPI_HOST_DEVICE HeatedCrossSectionMultiGroup( );
        LUPI_HOST HeatedCrossSectionMultiGroup( LUPI::StatusMessageReporting &a_smr, GIDI::ProtareSingle const &a_protare, SetupInfo &a_setupInfo, 
                Transporting::MC const &a_settings, GIDI::Styles::TemperatureInfo const &a_temperatureInfo,
                GIDI::Transporting::Particles const &a_particles, std::vector<GIDI::Reaction const *> const &a_reactions, std::string const &a_label,
                bool a_zeroReactions, GIDI::ExcludeReactionsSet const &a_reactionsToExclude );
        LUPI_HOST_DEVICE ~HeatedCrossSectionMultiGroup( );

        LUPI_HOST_DEVICE HeatedReactionCrossSectionMultiGroup *operator[]( std::size_t a_index ) const { return( m_reactionCrossSections[a_index] ); }
                                                                                /**< Returns the HeatedReactionCrossSectionMultiGroup for the reaction at index *a_index *a_index*. */
        LUPI_HOST_DEVICE int numberOfReactions( ) const { return( (int) m_reactionCrossSections.size( ) ); }
                                                                                /**< Returns the number of reactions stored in *this*. */

        LUPI_HOST_DEVICE int thresholdOffset(                              int a_index ) const { return( m_reactionCrossSections[a_index]->offset( ) ); }
                                                                                /**< Returns the offset for the cross section for the reaction with index *a_index*. */
        LUPI_HOST_DEVICE double threshold(                                 int a_index ) const { return( m_reactionCrossSections[a_index]->threshold( ) ); }

        LUPI_HOST_DEVICE Vector<double> &totalCrossSection( ) { return( m_totalCrossSection ); }    /**< Returns a reference to member *m_totalCrossSection*. */
        LUPI_HOST_DEVICE double crossSection(                              int a_hashIndex, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE double augmentedCrossSection(                     int a_hashIndex ) const { return( m_augmentedCrossSection[a_hashIndex] ); }
                                                                                /**< Returns the value of the of the augmented cross section the reaction at index *a_index*. */
        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, int a_hashIndex, bool a_sampling = false ) const {
                return( m_reactionCrossSections[a_reactionIndex]->crossSection( a_hashIndex, a_sampling ) ); }
                /**< Returns the reaction's cross section for the reaction at index *a_reactionIndex* and multi-group index *a_hashIndex*. */

        LUPI_HOST_DEVICE double depositionEnergy(   int a_hashIndex ) const { return( m_depositionEnergy[a_hashIndex] ); }
        LUPI_HOST_DEVICE double depositionMomentum( int a_hashIndex ) const { return( m_depositionMomentum[a_hashIndex] ); }
        LUPI_HOST_DEVICE double productionEnergy(   int a_hashIndex ) const { return( m_productionEnergy[a_hashIndex] ); }
        LUPI_HOST_DEVICE double gain(               int a_hashIndex, int a_particleIndex ) const ;
        LUPI_HOST_DEVICE double gainViaIntid(       int a_hashIndex, int a_particleIntid ) const ;

        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST void write( FILE *a_file ) const ;
};

/*
============================================================
============== HeatedCrossSectionsMultiGroup ==============
============================================================
*/
class HeatedCrossSectionsMultiGroup {

    private:
        Vector<double> m_temperatures;
        Vector<double> m_thresholds;
        Vector<int> m_multiGroupThresholdIndex;                         /**< This is the group where threshold starts, -1 otherwise. */
        Vector<double> m_projectileMultiGroupBoundariesCollapsed;
        Vector<HeatedCrossSectionMultiGroup *> m_heatedCrossSections;

    public:
        LUPI_HOST_DEVICE HeatedCrossSectionsMultiGroup( );
        LUPI_HOST_DEVICE ~HeatedCrossSectionsMultiGroup( );

        LUPI_HOST_DEVICE double minimumEnergy( ) const { return( m_projectileMultiGroupBoundariesCollapsed[0] ); }
        LUPI_HOST_DEVICE double maximumEnergy( ) const { return( m_projectileMultiGroupBoundariesCollapsed.back( ) ); }
        LUPI_HOST_DEVICE Vector<double> const &temperatures( ) const { return( m_temperatures ); }   /**< Returns the value of the **m_temperatures**. */

        LUPI_HOST void update( LUPI::StatusMessageReporting &a_smr, GIDI::ProtareSingle const &a_protare, SetupInfo &a_setupInfo, Transporting::MC const &a_settings, GIDI::Transporting::Particles const &a_particles, 
                GIDI::Styles::TemperatureInfos const &a_temperatureInfos, std::vector<GIDI::Reaction const *> const &a_reactions, 
                std::vector<GIDI::Reaction const *> const &a_orphanProducts, bool a_zeroReactions, GIDI::ExcludeReactionsSet const &a_reactionsToExclude );

        LUPI_HOST_DEVICE int multiGroupThresholdIndex( std::size_t a_index ) const { return( m_multiGroupThresholdIndex[a_index] ); }
                                                                                                    /**< Returns the threshold for the reaction at index *a_index*. */
        LUPI_HOST_DEVICE Vector<double> const &projectileMultiGroupBoundariesCollapsed( ) const { return( m_projectileMultiGroupBoundariesCollapsed ); }
                                                                                                    /**< Returns the value of the **m_projectileMultiGroupBoundariesCollapsed**. */
        LUPI_HOST_DEVICE Vector<HeatedCrossSectionMultiGroup *> const &heatedCrossSections( ) const { return( m_heatedCrossSections ); }

        LUPI_HOST_DEVICE double threshold( std::size_t a_index ) const { return( m_thresholds[a_index] ); }     /**< Returns the threshold for the reaction at index *a_index*. */

        LUPI_HOST_DEVICE double crossSection(                              int a_hashIndex, double a_temperature, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE void crossSectionVector( double a_temperature, double a_userFactor, std::size_t a_numberAllocated, 
                double *a_crossSectionVector ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, int a_hashIndex, double a_temperature, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, double a_temperature, double a_energy_in ) const ;
        template <typename RNG>
        inline LUPI_HOST_DEVICE int sampleReaction(                               int a_hashIndex, double a_temperature, double a_energy_in, double a_crossSection, RNG &&rng) const;

        LUPI_HOST_DEVICE double depositionEnergy(   int a_hashIndex, double a_temperature ) const ;
        LUPI_HOST_DEVICE double depositionMomentum( int a_hashIndex, double a_temperature ) const ;
        LUPI_HOST_DEVICE double productionEnergy(   int a_hashIndex, double a_temperature ) const ;
        LUPI_HOST_DEVICE double gain(               int a_hashIndex, double a_temperature, int a_particleIndex ) const ;
        LUPI_HOST_DEVICE double gainViaIntid(       int a_hashIndex, double a_temperature, int a_particleIntid ) const ;

        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST void write( FILE *a_file, int a_temperatureIndex ) const ;
        LUPI_HOST void print( ) const ;
};

/*
============================================================
================== NuclideGammaBranchInfo ==================
============================================================
*/
class NuclideGammaBranchInfo {

    private:
        double m_probability;                               /**< The probability that the level decays to state *m_residualStateIndex*. */
        double m_photonEmissionProbability;                 /**< The conditional probability the the decay emitted a photon. */
        double m_gammaEnergy;                               /**< The energy of the emitted photon. */
        int m_residualStateIndex;                           /**< The state the residual is left in after photon decay. */
        bool m_residualStateKindIsContinuum;                /**< True if the kind of the residual state (i.e., nuclide) is 'continuum' and false otherwise. */

    public:
        LUPI_HOST_DEVICE NuclideGammaBranchInfo( );
        LUPI_HOST NuclideGammaBranchInfo( PoPI::NuclideGammaBranchInfo const &a_nuclideGammaBranchInfo, 
                std::map<std::string, int> &a_stateNamesToIndices, bool a_makePhotonEmissionProbabilitiesOne );

        LUPI_HOST_DEVICE double probability( ) const { return( m_probability ); }                                  /**< Returns the value of the **m_probability**. */
        LUPI_HOST_DEVICE double photonEmissionProbability( ) const { return( m_photonEmissionProbability ); }      /**< Returns the value of the **m_photonEmissionProbability**. */
        LUPI_HOST_DEVICE double gammaEnergy( ) const { return( m_gammaEnergy ); }                                  /**< Returns the value of the **m_gammaEnergy**. */
        LUPI_HOST_DEVICE int residualStateIndex( ) const { return( m_residualStateIndex ); }                       /**< Returns the value of the **m_residualStateIndex**. */
        LUPI_HOST_DEVICE bool residualStateKindIsContinuum( ) const { return( m_residualStateKindIsContinuum ); }  /**< Returns the value of the **m_residualStateKindIsContinuum. **. */

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST void print( ProtareSingle const *a_protareSingle, std::string const &a_indent, std::string const &a_iFormat,
                std::string const &a_energyFormat, std::string const &a_dFormat ) const ;
};

/*
============================================================
============== NuclideGammaBranchStateInfo =================
============================================================
*/
class NuclideGammaBranchStateInfo {

    private:
        char m_state[16];                                   /**< The GNDS PoPs id for the nuclide. */
        int m_intid;                                        /**< The GNDS PoPs intid for the nuclide. */
        double m_nuclearLevelEnergy;                        /**< The nuclear level excitation energy of the level (state). */
        double m_nuclearLevelEnergyWidth;                   /**< This is 0.0 except for GRIN realized continuum levels where this is the energy width from this level to the next higher level. */
        double m_multiplicity;                              /**< The average multiplicity of photons emitted including the emission from sub-levels. */
        double m_averageGammaEnergy;                        /**< The average energy of photons emitted including the emission from sub-levels. */
        Vector<int> m_branchIndices;                        /**< The list of indices into the ProtareSingle.m_branches member that this level decays to. */

    public:
        LUPI_HOST_DEVICE NuclideGammaBranchStateInfo( );
        LUPI_HOST NuclideGammaBranchStateInfo( PoPI::NuclideGammaBranchStateInfo const &a_nuclideGammaBranchingInfo, 
                std::vector<NuclideGammaBranchInfo *> &a_nuclideGammaBranchInfos, 
                std::map<std::string, int> &a_stateNamesToIndices, bool a_makePhotonEmissionProbabilitiesOne,
                bool a_ignoreNuclearLevelEnergy );

        LUPI_HOST_DEVICE char const *state( ) const { return( m_state ); }                                    /**< Returns a pointer to the **m_state** member. */
        LUPI_HOST_DEVICE int intid( ) const { return( m_intid ); }                                            /**< Returns a pointer to the **m_intid** member. */
        LUPI_HOST_DEVICE double nuclearLevelEnergy( ) const { return( m_nuclearLevelEnergy ); }               /**< Returns the value of the **m_nuclearLevelEnergy** member. */
        LUPI_HOST_DEVICE double nuclearLevelEnergyWidth( ) const { return( m_nuclearLevelEnergyWidth ); }
                                                                                /**< Returns the value of the *m_nuclearLevelEnergyWidth* member. */
        LUPI_HOST_DEVICE double multiplicity( ) const { return( m_multiplicity ); }                           /**< Returns the value of the **m_multiplicity** member. */
        LUPI_HOST_DEVICE double averageGammaEnergy( ) const { return( m_averageGammaEnergy ); }               /**< Returns the value of the **m_averageGammaEnergy** member. */
        LUPI_HOST_DEVICE Vector<int> const &branchIndices( ) const { return( m_branchIndices ); }             /**< Returns the value of the **m_branchIndices** member. */

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST void print( ProtareSingle const *a_protareSingle, std::string const &a_indent, std::string const &a_iFormat,
                std::string const &a_energyFormat, std::string const &a_dFormat ) const ;
};

/*
============================================================
=============== GRIN_levelsAndProbabilities ================
============================================================
*/

class GRIN_levelsAndProbabilities {

    public:
        Vector<int> m_levels;                       /**< The list of nuclide indices for the nuclides in *m_state* as stored in member ProtareSingle::m_nuclideGammaBranchStateInfos. */
        Vector<double> m_summedProbabilities;       /**< The running sum of the probability for choosing a state from m_states. */
        Vector<bool> m_isModelledLevel;             /**< The entry for each item in *m_levels* which is true if the level is a modelled level and false otherwise. */

    public:
        LUPI_HOST_DEVICE  GRIN_levelsAndProbabilities( );
        LUPI_HOST         GRIN_levelsAndProbabilities( SetupInfo &a_setupInfo, PoPI::Database const &a_pops, 
                                GIDI::Table::Table const &a_table, bool a_normalize );
        LUPI_HOST_DEVICE ~GRIN_levelsAndProbabilities( );

        LUPI_HOST void set( std::vector<int> const &a_levels, std::vector<double> const &a_probabilities );

        template <typename RNG>
        inline LUPI_HOST_DEVICE int sampleInelasticLevel( double a_energy, RNG && a_rng );
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
================= GRIN_inelasticForEnergy ==================
============================================================
*/

class GRIN_inelasticForEnergy {

    private:
        Vector<int> m_indices;
        Vector<double> m_thresholds;
        GRIN_levelsAndProbabilities m_levelsAndProbabilities;

    public:
        LUPI_HOST_DEVICE GRIN_inelasticForEnergy( );
        LUPI_HOST GRIN_inelasticForEnergy( SetupInfo &a_setupInfo, double a_projectileMass, double a_targetMass,
                PoPI::Database const &a_pops, GIDI::GRIN::InelasticIncidentEnergy const *inelasticIncidentEnergy );
        LUPI_HOST_DEVICE ~GRIN_inelasticForEnergy( );

        LUPI_HOST_DEVICE int sampleLevelIndex( double a_projectileEnergy, double a_random ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
====================== GRIN_inelastic ======================
============================================================
*/

class GRIN_inelastic {

    private:
        int m_neutronIndex;
        int m_neutronUserParticleIndex;
        double m_neutronMass;

        int m_targetIntid;
        int m_targetIndex;
        int m_targetUserParticleIndex;
        double m_targetMass;

        Vector<double> m_energies;
        Vector<GRIN_inelasticForEnergy *> m_inelasticForEnergy;

    public:
        LUPI_HOST_DEVICE  GRIN_inelastic( );
        LUPI_HOST         GRIN_inelastic( SetupInfo &a_setupInfo, GIDI::GRIN::GRIN_continuumGammas const &GRIN_continuumGammas );
        LUPI_HOST_DEVICE ~GRIN_inelastic( );

        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );

        template <typename RNG, typename PUSHBACK>
        inline LUPI_HOST_DEVICE bool sampleProducts( ProtareSingle const *a_protare, double a_projectileEnergy, Sampling::Input &a_input, 
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
================= GRIN_captureToCompound ================
============================================================
*/

class GRIN_captureToCompound {

    private:
        int m_index;                                        /**< This is the index into ProtareSingle.m_nuclideGammaBranchStateInfos of the compound level forms by the capture. */
        GRIN_levelsAndProbabilities m_continuumIndices;     /**< This is the list of the levels the compound can decay to minus the known levels. */

    public:
        LUPI_HOST_DEVICE GRIN_captureToCompound( );
        LUPI_HOST GRIN_captureToCompound( SetupInfo &a_setupInfo, PoPI::Database const &a_pops, std::string a_compoundId );
        LUPI_HOST_DEVICE ~GRIN_captureToCompound( );

        LUPI_HOST_DEVICE int index( ) const { return( m_index ); }
        template <typename RNG>
        inline LUPI_HOST_DEVICE int sampleCaptureLevel( ProtareSingle const *a_protare, double a_energy, RNG && a_rng, bool a_checkEnergy ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
=============== GRIN_captureLevelProbability ===============
============================================================
*/

/* 
* The m_capturePrimaryToContinua are the modelled primary capture nuclides needed when one of the nuclide in 
* *m_levelsAndProbabilities* is not selected. The nuclide in *m_capturePrimaryToContinua* whose nuclear level 
* energy as found in *m_nuclearLevelEnergies* is closest but below the * neutron separation energy plus the 
* kinetic energy in the com frame is to be used when a nuclide in *m_levelsAndProbabilities* is not selected.
* the ints in *m_capturePrimaryToContinua* are indicies into ProtareSingle.m_nuclideGammaBranchStateInfos.
*/

class GRIN_captureLevelProbability {

    private:
        GRIN_levelsAndProbabilities m_knownLevelsAndProbabilities;      /**< These are the known primary capture nuclide as probabilites. The probabilites do not sum to one since not all are know. */
        Vector<GRIN_captureToCompound *> m_captureToCompounds;          /**< This is a list of nuclides that have the same spin/parity as the input channel and the needed nuclear excitation level. */

    public:
        LUPI_HOST_DEVICE  GRIN_captureLevelProbability( );
        LUPI_HOST         GRIN_captureLevelProbability( SetupInfo &a_setupInfo, PoPI::Database const &a_pops, 
                GIDI::GRIN::CaptureLevelProbability const *a_captureLevelProbability );
        LUPI_HOST_DEVICE ~GRIN_captureLevelProbability( );

        template <typename RNG>
        inline LUPI_HOST_DEVICE int sampleCaptureLevel( ProtareSingle const *a_protare, double a_energy, RNG && a_rng );
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
====================== GRIN_capture ========================
============================================================
*/

class GRIN_capture {

    private:
        double m_captureNeutronSeparationEnergy;                                /**< For capture, the neutron separation energy as needed to emit primary gammas.  */
        Vector<double> m_summedProbabilities;                                   /**< The running sum of the probabilites for the *m_captureLevelProbabilities* member data. */
        Vector<GRIN_captureLevelProbability *> m_captureLevelProbabilities;     /**< The list of capture levels with probabilites. */
        int m_residualIntid;                                                    /**< The intid of the heavy residual particle. */
        int m_residualIndex;                                                    /**< The PoPI index of the heavy residual particle. */
        int m_residualUserIndex;                                                /**< The user index of the heavy residual particle. */
        double m_residualMass;                                                  /**< The mass if the heavy residual particle. */

    public:
        LUPI_HOST_DEVICE  GRIN_capture( );
        LUPI_HOST         GRIN_capture( SetupInfo &a_setupInfo, GIDI::GRIN::GRIN_continuumGammas const &GRIN_continuumGammas );
        LUPI_HOST_DEVICE ~GRIN_capture( );

        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );
        template <typename RNG, typename PUSHBACK>
        inline LUPI_HOST_DEVICE bool sampleProducts( ProtareSingle const *a_protare, double a_projectileEnergy, Sampling::Input &a_input, 
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================= Product ==========================
============================================================
*/
class Product {

    private:
        String m_ID;
        int m_intid;
        int m_index;
        int m_userParticleIndex;
        String m_label;
        bool m_isCompleteParticle;
        double m_mass;
        double m_excitationEnergy;
        TwoBodyOrder m_twoBodyOrder;
        int m_initialStateIndex;                                /**< If the product has branching photons, then this is the state index to start the branching. */
        Functions::Function1d *m_multiplicity;
        Distributions::Distribution *m_distribution;
// still need *m_averageEnergy *m_averageMomentum;

        OutputChannel *m_outputChannel;

    public:
        LUPI_HOST_DEVICE Product( );
        LUPI_HOST Product( GIDI::Product const *a_product, SetupInfo &a_setupInfo, Transporting::MC const &a_settings, GIDI::Transporting::Particles const &a_particles,
                bool a_isFission );
        LUPI_HOST Product( PoPI::Database const &a_pop, std::string const &a_ID, std::string const &a_label );
        LUPI_HOST_DEVICE ~Product( );

        LUPI_HOST String const &ID( ) const { return( m_ID ); }                                 /**< Returns a const reference to the *m_ID* member. */
        LUPI_HOST_DEVICE int intid( ) const { return( m_intid); }                               /**< Returns the value of the *m_intid* member. */
        LUPI_HOST_DEVICE int index( ) const { return( m_index); }                               /**< Returns the value of the *m_index* member. */
        LUPI_HOST_DEVICE int userParticleIndex( ) const { return( m_userParticleIndex ); }      /**< Returns the value of the **m_userParticleIndex**. */
        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );
        LUPI_HOST void setModelDBRC_data( Sampling::Upscatter::ModelDBRC_data *a_modelDBRC_data );
        LUPI_HOST_DEVICE String label( ) const { return( m_label ); }                           /**< Returns the value of the **m_label**. */
        LUPI_HOST_DEVICE bool isCompleteParticle( ) const { return( m_isCompleteParticle ); }   /**< Returns the value of the **m_isCompleteParticle**. */
        LUPI_HOST_DEVICE double mass( ) const { return( m_mass ); }                             /**< Returns the value of the **m_mass**. */
        LUPI_HOST_DEVICE double excitationEnergy( ) const { return( m_excitationEnergy ); }     /**< Returns the value of the **m_excitationEnergy**. */
        LUPI_HOST_DEVICE TwoBodyOrder twoBodyOrder( ) const { return( m_twoBodyOrder ); }       /**< Returns the value of the **m_twoBodyOrder**. */
        LUPI_HOST_DEVICE double finalQ( double a_x1 ) const ;
        LUPI_HOST_DEVICE bool hasFission( ) const ;

// FIXME (1) see FIXME (1) in MC class.
        LUPI_HOST_DEVICE Functions::Function1d const *multiplicity( ) const { return( m_multiplicity ); }      /**< Returns the value of the **m_multiplicity**. */
        LUPI_HOST void setMultiplicity( Functions::Function1d *a_multiplicity ) { m_multiplicity = a_multiplicity; }
        LUPI_HOST_DEVICE double productAverageMultiplicity(         int a_index, double a_projectileEnergy ) const ;
        LUPI_HOST_DEVICE double productAverageMultiplicityViaIntid( int a_intid, double a_projectileEnergy ) const ;
// FIXME (1) see FIXME (1) in MC class.
        LUPI_HOST_DEVICE Distributions::Distribution const *distribution( ) const { return( m_distribution ); }     /**< Returns the value of the **m_distribution**. */
        LUPI_HOST_DEVICE Distributions::Distribution *distribution( ) { return( m_distribution ); }                 /**< Returns the value of the **m_distribution**. */
        LUPI_HOST void distribution( Distributions::Distribution *a_distribution ) { m_distribution = a_distribution; }
// FIXME (1) see FIXME (1) in MC class.
        LUPI_HOST_DEVICE OutputChannel *outputChannel( ) { return( m_outputChannel ); }                  /**< Returns the value of the **m_outputChannel**. */

        template <typename RNG, typename PUSHBACK>
        inline LUPI_HOST_DEVICE void sampleProducts( ProtareSingle const *a_protare, double a_projectileEnergy, Sampling::Input &a_input,
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const ;
        template <typename RNG, typename PUSHBACK>
        inline LUPI_HOST_DEVICE void sampleFinalState( ProtareSingle const *a_protare, double a_projectileEnergy, Sampling::Input &a_input,
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const ;
        template <typename RNG>
        inline LUPI_HOST_DEVICE void angleBiasing( Reaction const *a_reaction, int a_pid, double a_temperature, double a_energy_in, double a_mu_lab, 
              double &a_probability, double &a_energy_out, RNG && a_rng, double &a_cumulative_weight ) const ;
        template <typename RNG>
        inline LUPI_HOST_DEVICE void angleBiasingViaIntid( Reaction const *a_reaction, int a_intid, double a_temperature, double a_energy_in, double a_mu_lab, 
                double &a_probability, double &a_energy_out, RNG && a_rng, double &a_cumulative_weight ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
====================== DelayedNeutron ======================
============================================================
*/
class DelayedNeutron {

    private:
        int m_delayedNeutronIndex;                  /**< If this is a delayed fission neutron, this is its index. */
        double m_rate;                              /**< The GNDS rate for the delayed neutron. */
        Product m_product;                          /**< The GNDS <**product**> node. */

    public:
        LUPI_HOST_DEVICE DelayedNeutron( );
        LUPI_HOST DelayedNeutron( int a_index, GIDI::DelayedNeutron const *a_delayedNeutron, SetupInfo &a_setupInfo, Transporting::MC const &a_settings, GIDI::Transporting::Particles const &a_particles );
        LUPI_HOST_DEVICE ~DelayedNeutron( );

        LUPI_HOST_DEVICE int delayedNeutronIndex( ) const { return( m_delayedNeutronIndex ); }
        LUPI_HOST_DEVICE double rate( ) const { return( m_rate ); }
        LUPI_HOST_DEVICE Product const &product( ) const { return( m_product ); }
        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======================= OutputChannel ======================
============================================================
*/
class OutputChannel {

    private:
        ChannelType m_channelType;
        int m_neutronIndex;                         /**< The index of the neutron in the use PoPs database. */
        bool m_isFission;
        bool m_hasFinalStatePhotons;                /**< If **true**, *this* channel has a photon with finalState attribute. */

        Functions::Function1d_d1 *m_Q;              /**< The Q-function for the output channel. Note, this is currently always the *evaluated* form even when running with multi-group data. */
        Vector<Product *> m_products;
        Functions::Function1d *m_totalDelayedNeutronMultiplicity;
        Vector<DelayedNeutron *> m_delayedNeutrons;

    public:
        LUPI_HOST_DEVICE OutputChannel( );
        LUPI_HOST OutputChannel( GIDI::OutputChannel const *a_outputChannel, SetupInfo &a_setupInfo, Transporting::MC const &a_settings, GIDI::Transporting::Particles const &a_particles );
        LUPI_HOST_DEVICE ~OutputChannel( );

        LUPI_HOST_DEVICE Product *operator[]( std::size_t a_index ) { return( m_products[a_index] ); }  /**< Returns a pointer to the product at index *a_index*. */

        LUPI_HOST_DEVICE bool isTwoBody( ) const { return( m_channelType == ChannelType::twoBody ); }         /**< Returns true if output channel is two-body and false otherwise. */
        LUPI_HOST_DEVICE double finalQ( double a_x1 ) const ;
        LUPI_HOST_DEVICE bool isFission( ) const { return( m_isFission ); }                      /**< Returns the value of the **m_isFission**. */
        LUPI_HOST_DEVICE bool hasFission( ) const ;
// FIXME (1) see FIXME (1) in MC class.
        LUPI_HOST_DEVICE Functions::Function1d_d1 *Q( ) { return( m_Q ); }                       /**< Returns the pointer of the **m_Q** member. */

        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );
        LUPI_HOST void setModelDBRC_data( Sampling::Upscatter::ModelDBRC_data *a_modelDBRC_data );
// FIXME (1) see FIXME (1) in MC class.
        LUPI_HOST_DEVICE Vector<Product *> const &products( ) const { return( m_products ); }    /**< Returns the value of the **m_products**. */

        Vector<DelayedNeutron *> delayedNeutrons( ) const { return( m_delayedNeutrons ); }
        LUPI_HOST_DEVICE DelayedNeutron const *delayedNeutron( int a_index ) const { return( m_delayedNeutrons[a_index] ); }

        LUPI_HOST void moveProductsEtAlToReaction( std::vector<Product *> &a_products, Functions::Function1d **a_totalDelayedNeutronMultiplicity, 
                std::vector<DelayedNeutron *> &a_delayedNeutrons, std::vector<Functions::Function1d_d1 *> &a_Qs );
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
        LUPI_HOST void addOrphanProductToProductList( std::vector<Product *> &a_associatedOrphanProducts ) const ;
        LUPI_HOST_DEVICE void addOrphanProductToProductList( Vector<Product *> &a_associatedOrphanProducts ) const ;
#endif

        LUPI_HOST_DEVICE double productAverageMultiplicity(         int a_index, double a_projectileEnergy ) const ;
        LUPI_HOST_DEVICE double productAverageMultiplicityViaIntid( int a_intid, double a_projectileEnergy ) const ;

template <typename RNG, typename PUSHBACK>
        inline LUPI_HOST_DEVICE void sampleProducts( ProtareSingle const *a_protare, double a_projectileEnergy, Sampling::Input &a_input,
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const;
        template <typename RNG>
        inline LUPI_HOST_DEVICE void angleBiasing( Reaction const *a_reaction, int a_pid, double a_temperature, double a_energy_in, double a_mu_lab, 
              double &a_probability, double &a_energy_out, RNG && a_rng, double &a_cumulative_weight ) const ;
        template <typename RNG>
        inline LUPI_HOST_DEVICE void angleBiasingViaIntid( Reaction const *a_reaction, int a_intid, double a_temperature, double a_energy_in, double a_mu_lab, 
                double &a_probability, double &a_energy_out, RNG && a_rng, double &a_cumulative_weight ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================= Reaction =========================
============================================================
*/
class Reaction {

    private:
        ProtareSingle *m_protareSingle;                     /**< The ProtareSingle this reaction resides in. */
        int m_reactionIndex;                                /**< The index of the reaction in the ProtareSingle. */
        int m_GIDI_reactionIndex;                           /**< The index of the reaction in the GIDI::ProtareSingle. */
        String m_label;                                     /**< The **GNDS** label for the reaction. */
        int m_ENDF_MT;                                      /**< The ENDF MT value for the reaction. */
        int m_ENDL_C;                                       /**< The ENDL C value for the reaction. */
        int m_ENDL_S;                                       /**< The ENDL S value for the reaction. */
        int m_initialStateIndex;                            /**< If not -1, then reaction contains a branching gamma data with this index for the data in m_nuclideGammaBranchStateInfos member of its **ProtareSingle** instance. */
        int m_neutronIndex;                                 /**< The index of the neutron in the use PoPs database. */
        bool m_hasFission;                                  /**< Is *true* if the reaction is a fission reaction and *false* otherwise. */
        double m_projectileMass;                            /**< The mass of the projectile. */
        double m_targetMass;                                /**< The mass of the target. */
        double m_crossSectionThreshold;                     /**< The threshold for the reaction. */
        double m_twoBodyThreshold;                          /**< This is the T_1 value needed to do two-body kinematics. */
        bool m_upscatterModelASupported;
        bool m_hasFinalStatePhotons;                        /**< If **true**, *this* reaction has a photon with finalState attribute. */
        int m_fissionResiduaIntid;                          /**< The intid of the special ENDL 99120 or 99125 fission residual. */
        int m_fissionResiduaIndex;                          /**< The index of the special ENDL 99120 or 99125 fission residual. */
        int m_fissionResiduaUserIndex;                      /**< The user index of the special ENDL 99120 or 99125 fission residual. */
        GIDI::Construction::FissionResiduals m_fissionResiduals;  /**< This member specifies what fission redisual products will be added to the list of products produced in a fission reaction. */
        double m_fissionResidualMass;                       /**< The mass of the special ENDL 99120 or 99125 fission residual. */
        Vector<double> m_upscatterModelACrossSection;       /**< The multi-group cross section to use for upscatter model A. */

        Vector<int> m_productIntids;                        /**< The list of all products *this* reaction can product by their intid. */
        Vector<int> m_productIndices;                       /**< The list of all products *this* reaction can product by their index. */
        Vector<int> m_userProductIndices;                   /**< The list of all products *this* reaction can product as user indices. */
        Vector<int> m_productMultiplicities;                /**< The list of all multiplicities for each product in *m_productIntids* . */
        Vector<int> m_productIntidsTransportable;           /**< The list of all transportabls products *this* reaction can product by their intid. */
        Vector<int> m_productIndicesTransportable;          /**< The list of all transportabls products *this* reaction can product by their index. */
        Vector<int> m_userProductIndicesTransportable;      /**< The list of all transportabls products *this* reaction can product as user indices. */

        Vector<Functions::Function1d_d1 *> m_Qs;            /**< A list of Q-functions that is used when the C macro MCGIDI_USE_OUTPUT is defined. */
        Vector<Product *> m_products;                       /**< A list of all transporting products directly or nested in **m_outputChannel** that is used instead of having **m_outputChannel** loop of all transporting products if the C macro MCGIDI_USE_OUTPUT_CHANNEL is not defined. */
        Functions::Function1d *m_totalDelayedNeutronMultiplicity;
        Vector<DelayedNeutron *> m_delayedNeutrons;         /**< A list of all delayedNeutrons that can be used instead of having m_outputChannel loop of all transporting products. For *m_products* for more details. */
                                                            /**< The total delayed neutron multiplicity used when the C macro MCGIDI_USE_OUTPUT is defined. */
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
        OutputChannel *m_outputChannel;                     /**< The output channel for this reaction. Only used if the C macro MCGIDI_USE_OUTPUT is defined. */
#endif
        Vector<int> m_associatedOrphanProductIndices;       /**< The indices in the Protare's m_orphanProducts member for the orphanProducts associated with this reaction. */
        Vector<Product *> m_associatedOrphanProducts;       /**< The list of products from the orphanProduct reaction. */ /* Do not delete entries as owned by orphanProduct reaction. */
// Still need m_availableEnergy and m_availableMomentum.

// GRIN specials: new non-GNDS GRIN stuff.
        bool m_GRIN_specialSampleProducts;                  /**< This will be true if sampling products with GRIN special inelastic or capture data. */
        double m_GRIN_inelasticThreshold;                   /**< For inelastic, the a_projectileEnergy must be greater than this value 
                                                                or use standard product sampling.  This is needed, for example, as the Fe56 MT 91 cross section 
                                                                starts below the MT 89 cross section. Ergo, below the threshold for the nuclear level energy 
                                                                for the first simulated Fe56 nuclear level. */
        double m_GRIN_maximumCaptureIncidentEnergy;         /**< For capture, the projectile energy must be less thans this value or use standard product sampling. */
        GRIN_inelastic *m_GRIN_inelastic;                   /**< A nullptr or a pointer to an instance of GIND_continuumInelatic (see below). */
        GRIN_capture *m_GRIN_capture;                       /**< A nullptr or a pointer to an instance of GIND_capture (see below).  */

    public:
        LUPI_HOST_DEVICE Reaction( );
        LUPI_HOST Reaction( GIDI::Reaction const &a_reaction, SetupInfo &a_setupInfo, Transporting::MC const &a_settings, GIDI::Transporting::Particles const &a_particles,
                GIDI::Styles::TemperatureInfos const &a_temperatureInfos );
        LUPI_HOST_DEVICE ~Reaction( );

        inline LUPI_HOST_DEVICE void updateProtareSingleInfo( ProtareSingle *a_protareSingle, int a_reactionIndex ) {
                m_protareSingle = a_protareSingle;
                m_reactionIndex = a_reactionIndex;
        }
        LUPI_HOST_DEVICE ProtareSingle const *protareSingle( ) const { return( m_protareSingle ); }   /**< Returns the value of the **m_protareSingle**. */
        LUPI_HOST_DEVICE int reactionIndex( ) const { return( m_reactionIndex ); }       /**< Returns the value of the **m_reactionIndex**. */
        LUPI_HOST_DEVICE int GIDI_reactionIndex( ) const { return( m_GIDI_reactionIndex ); }       /**< Returns the value of the **m_GIDI_reactionIndex** member. */
        LUPI_HOST_DEVICE String const &label( ) const { return( m_label ); }             /**< Returns the value of the **m_label**. */
        LUPI_HOST_DEVICE int ENDF_MT( ) const { return( m_ENDF_MT ); }                   /**< Returns the value of the **m_ENDF_MT**. */
        LUPI_HOST_DEVICE int ENDL_C( ) const { return( m_ENDL_C ); }                     /**< Returns the value of the **m_ENDL_C**. */
        LUPI_HOST_DEVICE int ENDL_S( ) const { return( m_ENDL_S ); }                     /**< Returns the value of the **m_ENDL_S**. */
        LUPI_HOST_DEVICE int initialStateIndex( ) const { return( m_initialStateIndex ); }    /**< Returns the value of the **m_initialStateIndex** member. */
        LUPI_HOST_DEVICE double finalQ( double a_energy ) const ;
        LUPI_HOST_DEVICE bool hasFission( ) const { return( m_hasFission ); }            /**< Returns the value of the **m_hasFission**. */
        LUPI_HOST_DEVICE double projectileMass( ) const { return( m_projectileMass ); }  /**< Returns the value of the **m_projectileMass**. */
        LUPI_HOST_DEVICE double targetMass( ) const { return( m_targetMass ); }          /**< Returns the value of the **m_targetMass**. */
        LUPI_HOST_DEVICE double crossSectionThreshold( ) const { return( m_crossSectionThreshold ); }    /**< Returns the value of the **m_crossSectionThreshold**. */
        LUPI_HOST_DEVICE double twoBodyThreshold( ) const { return( m_twoBodyThreshold ); }              /**< Returns the value of the *m_twoBodyThreshold* member. */
        LUPI_HOST_DEVICE double crossSection( URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double crossSection( URR_protareInfos const &a_URR_protareInfos, double a_temperature, double a_energy ) const ;
        LUPI_HOST GIDI::Functions::XYs1d crossSectionAsGIDI_XYs1d( double a_temperature ) const ;

        LUPI_HOST_DEVICE Vector<int> const &productIntids( ) const { return( m_productIntids ); }
        LUPI_HOST_DEVICE Vector<int> const &productIndices( ) const { return( m_productIndices ); }            /**< Returns a const reference to the *m_productIntids* member. */
        LUPI_HOST_DEVICE Vector<int> const &userProductIndices( ) const { return( m_userProductIndices ); }    /**< Returns a const reference to the *m_productIndices* member. */
        LUPI_HOST_DEVICE MCGIDI_VectorSizeType numberOfProducts( ) const { return( m_products.size( ) ); }     /**< Returns the number of products in the **m_products** member. */
        LUPI_HOST_DEVICE Product const *product( int a_index ) const { return( m_products[a_index] ); }
        LUPI_HOST_DEVICE int productMultiplicity(         int a_index ) const ;
        LUPI_HOST_DEVICE int productMultiplicityViaIntid( int a_intid ) const ;
        LUPI_HOST_DEVICE int productMultiplicities( int a_index ) const {
                LUPI::deprecatedFunction( "MCGIDI::Reaction::productMultiplicities", "MCGIDI::Reaction::productMultiplicity", "" );
                return( productMultiplicity( a_index ) ); }                                 /**< This method is deprecated. Please use **productMultiplicity** instead. */
        LUPI_HOST_DEVICE double productAverageMultiplicity(         int a_index, double a_projectileEnergy ) const ;
        LUPI_HOST_DEVICE double productAverageMultiplicityViaIntid( int a_intid, double a_projectileEnergy ) const ;
        LUPI_HOST_DEVICE Vector<int> const &productIntidsTransportable( ) const { return( m_productIntidsTransportable ); }
                                                                                            /**< Returns a const reference to the *m_productIntidsTransportable* member. */
        LUPI_HOST_DEVICE Vector<int> const &productIndicesTransportable( ) const { return( m_productIndicesTransportable ); }
                                                                                            /**< Returns a const reference to the *m_productIndicesTransportable* member. */
        LUPI_HOST_DEVICE Vector<int> const &userProductIndicesTransportable( ) const { return( m_userProductIndicesTransportable ); }

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
        LUPI_HOST_DEVICE OutputChannel const *outputChannel( ) const { return( m_outputChannel ); }              /**< Returns the value of the **m_outputChannel**. */
#endif
        LUPI_HOST_DEVICE Vector<int> associatedOrphanProductIndices( ) const { return( m_associatedOrphanProductIndices ); } /**< Returns the value of the **m_associatedOrphanProductIndicex** member. */
        LUPI_HOST void addOrphanProductToProductList( std::vector<Product *> &a_associatedOrphanProducts ) const ;
        LUPI_HOST_DEVICE void addOrphanProductToProductList( Vector<Product *> &a_associatedOrphanProducts ) const ;
        LUPI_HOST_DEVICE void addOrphanProductToProductList( Vector<Reaction *> &a_orphanProducts ) ;
        LUPI_HOST void setOrphanProductData( std::vector<int> const &a_associatedOrphanProductIndcies,
                std::vector<Product *> const &a_associatedOrphanProducts );
        LUPI_HOST_DEVICE bool upscatterModelASupported( ) const { return( m_upscatterModelASupported ); }
        LUPI_HOST_DEVICE Vector<double> const &upscatterModelACrossSection( ) const { return( m_upscatterModelACrossSection ); } 
                                                                                                            /**< Returns the value of the **m_upscatterModelACrossSection**. */

        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );
        LUPI_HOST void setModelDBRC_data( Sampling::Upscatter::ModelDBRC_data *a_modelDBRC_data );

        template <typename RNG, typename PUSHBACK>
        inline LUPI_HOST_DEVICE void sampleProducts( Protare const *a_protare, double a_projectileEnergy, Sampling::Input &a_input, 
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products, bool a_checkOrphanProducts = true ) const ;
        template <typename RNG, typename PUSHBACK>
        inline LUPI_HOST_DEVICE static void sampleNullProducts( Protare const &a_protare, double a_projectileEnergy, Sampling::Input &a_input, 
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products );
        template <typename RNG>
        inline LUPI_HOST_DEVICE double angleBiasing( int a_pid, double a_temperature, double a_energy_in, double a_mu_lab, double &a_energy_out, 
                RNG && a_rng, double *a_cumulative_weight = nullptr, bool a_checkOrphanProducts = true ) const ;
        template <typename RNG>
        inline LUPI_HOST_DEVICE double angleBiasingViaIntid( int a_intid, double a_temperature, double a_energy_in, double a_mu_lab, double &a_energy_out, 
                RNG && a_rng, double *a_cumulative_weight = nullptr, bool a_checkOrphanProducts = true ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================== Protare =========================
============================================================
*/
class Protare {

    private:
        ProtareType m_protareType;                          /**< The type of protare *this* is. */

        String m_projectileID;                              /**< The PoPs id of the projectile. */
        int m_projectileIntid;                              /**< The PoPs intid of the projectile. */
        int m_projectileIndex;                              /**< The PoPs database index of the projectile. */
        int m_projectileUserIndex;                          /**< The projectile's index as specified by the user. */
        double m_projectileMass;                            /**< The mass of the projectile. */
        double m_projectileExcitationEnergy;                /**< The nuclear excitation of the projectile. */

        String m_targetID;                                  /**< The PoPs intid of the target. */
        int m_targetIntid;                                  /**< The PoPs index of the target. */
        int m_targetIndex;                                  /**< The PoPs database index of the target. */
        int m_targetUserIndex;                              /**< The target's index as specified by the user. */
        double m_targetMass;                                /**< The mass of the target. */
        double m_targetExcitationEnergy;                    /**< The nuclear excitation of the target. */

        int m_neutronIndex;                                 /**< The neutron particle index from the user's pops database. */
        int m_userNeutronIndex;                             /**< The neutron particle index defined by the user. */
        int m_photonIndex;                                  /**< The photon particle index from the user's pops database. */
        int m_userPhotonIndex;                              /**< The photon particle index defined by the user. */

        String m_evaluation;                                /**< The evaluation string for the Protare. */
        GIDI::Frame m_projectileFrame;                      /**< The frame the projectile data are given in. */

        Vector<int> m_productIntids;                        /**< The list of all products *this* protare can product by their intid. */
        Vector<int> m_productIndices;                       /**< The list of all products *this* reaction can product by their index. */
        Vector<int> m_userProductIndices;                   /**< The list of all products *this* reaction can product as user indices. */
        Vector<int> m_productIntidsTransportable;           /**< The list of all transportabls products *this* protare can product by their intid. */
        Vector<int> m_productIndicesTransportable;          /**< The list of all transportabls products *this* reaction can product by their index. */
        Vector<int> m_userProductIndicesTransportable;      /**< The list of all transportabls products *this* reaction can product as user indices. */

        bool m_isTNSL_ProtareSingle;                        /**< If *this* is a ProtareSingle instance with TNSL data *true* and otherwise *false*. */

        LUPI_HOST void productIntidsAndIndices( std::set<int> const &a_intids, std::set<int> const &a_transportableIntids,
                std::set<int> const &a_indices, std::set<int> const &a_transportableIndices );

    public:
        LUPI_HOST_DEVICE Protare( ProtareType a_protareType );
        LUPI_HOST Protare( ProtareType a_protareType, GIDI::Protare const &a_protare, Transporting::MC const &a_settings, PoPI::Database const &a_pops );
        virtual LUPI_HOST_DEVICE ~Protare( );

        LUPI_HOST_DEVICE ProtareType protareType( ) const { return( m_protareType ); }                               /**< Returns the value of the **m_protareType** member. */    

        LUPI_HOST_DEVICE String const &projectileID( ) const { return( m_projectileID ); }                       /**< Returns the value of the **m_projectileID** member. */
        LUPI_HOST_DEVICE int projectileIntid( ) const { return( m_projectileIntid ); }                           /**< Returns the value of the **m_projectileIntid** member. */
        LUPI_HOST_DEVICE int projectileIndex( ) const { return( m_projectileIndex ); }                           /**< Returns the value of the **m_projectileIndex** member. */
        LUPI_HOST_DEVICE int projectileUserIndex( ) const { return( m_projectileUserIndex ); }                   /**< Returns the value of the **m_projectileUserIndex** member. */
        LUPI_HOST_DEVICE double projectileMass( ) const { return( m_projectileMass ); }                          /**< Returns the value of the **m_projectileMass** member. */
        LUPI_HOST_DEVICE double projectileExcitationEnergy( ) const { return( m_projectileExcitationEnergy ); }  /**< Returns the value of the **m_projectileExcitationEnergy** member. */

        LUPI_HOST_DEVICE String const &targetID( ) const { return( m_targetID ); }                               /**< Returns the value of the **m_targetID** member. */
        LUPI_HOST_DEVICE int targetIntid( ) const { return( m_targetIntid ); }                                   /**< Returns the value of the **m_targetIntid** member. */
        LUPI_HOST_DEVICE int targetIndex( ) const { return( m_targetIndex ); }                                   /**< Returns the value of the **m_targetIndex** member. */
        LUPI_HOST_DEVICE int targetUserIndex( ) const { return( m_targetUserIndex ); }                           /**< Returns the value of the **m_targetUserIndex** member. */
        LUPI_HOST_DEVICE double targetMass( ) const { return( m_targetMass ); }                                  /**< Returns the value of the **m_targetMass** member. */
        LUPI_HOST_DEVICE double targetExcitationEnergy( ) const { return( m_targetExcitationEnergy ); }          /**< Returns the value of the **m_targetExcitationEnergy** member. */

        LUPI_HOST_DEVICE int photonIndex( ) const { return( m_photonIndex ); }                                   /**< Returns the value of the **m_photonIndex** member. */
        LUPI_HOST_DEVICE int userPhotonIndex( ) const { return( m_userPhotonIndex ); }                           /**< Returns the value of the **m_userPhotonIndex** member. */
        LUPI_HOST_DEVICE String evaluation( ) const { return( m_evaluation ); }                                  /**< Returns the value of the **m_evaluation** member. */
        LUPI_HOST GIDI::Frame projectileFrame( ) const { return( m_projectileFrame ); }                          /**< Returns the value of the **m_projectileFrame** member. */

        LUPI_HOST Vector<int> const &productIntids( bool a_transportablesOnly ) const ;
        LUPI_HOST Vector<int> const &productIndices( bool a_transportablesOnly ) const ;
        LUPI_HOST Vector<int> const &userProductIndices( bool a_transportablesOnly ) const ;
        LUPI_HOST void setUserParticleIndex( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid( int a_particleIntid, int a_userParticleIndex );

        LUPI_HOST_DEVICE bool isTNSL_ProtareSingle( ) const { return( m_isTNSL_ProtareSingle ); }                /**< Returns the value of the **m_isTNSL_ProtareSingle** member. */
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE std::size_t numberOfProtares( ) const MCGIDI_TRUE_VIRTUAL;                            /**< Returns the number of protares contained in *this*. */
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE ProtareSingle const *protare( std::size_t a_index ) const MCGIDI_TRUE_VIRTUAL;        /**< Returns the **a_index** - 1 Protare contained in *this*. */
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE ProtareSingle       *protare( std::size_t a_index )       MCGIDI_TRUE_VIRTUAL;        /**< Returns the **a_index** - 1 Protare contained in *this*. */
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE ProtareSingle const *protareWithReaction( int a_index ) const MCGIDI_TRUE_VIRTUAL;              /**< Returns the *ProtareSingle* that contains the (*a_index* - 1) reaction. */

        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double minimumEnergy( ) const MCGIDI_TRUE_VIRTUAL;                                              /**< Returns the minimum cross section domain. */
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double maximumEnergy( ) const MCGIDI_TRUE_VIRTUAL ;                                             /**< Returns the maximum cross section domain. */
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE Vector<double> temperatures( std::size_t a_index = 0 ) const MCGIDI_TRUE_VIRTUAL ;    /**< Returns the list of temperatures for the requested ProtareSingle. */

        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST Vector<double> const &projectileMultiGroupBoundaries( ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST Vector<double> const &projectileMultiGroupBoundariesCollapsed( ) const MCGIDI_TRUE_VIRTUAL;

        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE std::size_t numberOfReactions( ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE Reaction const *reaction( int a_index ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE std::size_t numberOfOrphanProducts( ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE Reaction const *orphanProduct( int a_index ) const MCGIDI_TRUE_VIRTUAL;

        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE bool hasFission( ) const MCGIDI_TRUE_VIRTUAL;

        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE bool hasIncoherentDoppler( ) const MCGIDI_TRUE_VIRTUAL;

        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE int URR_index( ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE bool hasURR_probabilityTables( ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double URR_domainMin( ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double URR_domainMax( ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE bool reactionHasURR_probabilityTables( int a_index ) const MCGIDI_TRUE_VIRTUAL ;

        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double threshold( std::size_t a_index ) const MCGIDI_TRUE_VIRTUAL;

        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double crossSection(                              URR_protareInfos const &a_URR_protareInfos,
                int a_hashIndex, double a_temperature, double a_energy, bool a_sampling = false ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE void crossSectionVector( double a_temperature, double a_userFactor, int a_numberAllocated,
                double *a_crossSectionVector ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_hashIndex,
                double a_temperature, double a_energy, bool a_sampling = false ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos,
                double a_temperature, double a_energy ) const MCGIDI_TRUE_VIRTUAL;
        template <typename RNG>
        inline MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE int sampleReaction(                               URR_protareInfos const &a_URR_protareInfos, int a_hashIndex,
                double a_temperature, double a_energy, double a_crossSection, RNG && a_rng) const MCGIDI_TRUE_VIRTUAL;

        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double depositionEnergy(   int a_hashIndex, double a_temperature, double a_energy ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double depositionMomentum( int a_hashIndex, double a_temperature, double a_energy ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double productionEnergy(   int a_hashIndex, double a_temperature, double a_energy ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double gain(               int a_hashIndex, double a_temperature, double a_energy, int a_particleIndex ) const MCGIDI_TRUE_VIRTUAL;
        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE double gainViaIntid(       int a_hashIndex, double a_temperature, double a_energy, int a_particleIntid ) const MCGIDI_TRUE_VIRTUAL;

        MCGIDI_VIRTUAL_FUNCTION LUPI_HOST_DEVICE Vector<double> const &upscatterModelAGroupVelocities( ) const MCGIDI_TRUE_VIRTUAL;

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST_DEVICE void serialize2( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST_DEVICE void serializeCommon( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST_DEVICE long sizeOf( ) const ;
        LUPI_HOST_DEVICE long memorySize( );
        LUPI_HOST_DEVICE void incrementMemorySize( long &a_totalMemory, long &a_sharedMemory );

        friend ProtareSingle;
        friend ProtareComposite;
        friend ProtareTNSL;
};

/*
============================================================
====================== ProtareSingle =======================
============================================================
*/
class ProtareSingle : public Protare {

    private:
        String m_interaction;                                                       /**< The protare's interaction string. */
        int m_URR_index;                                                            /**< The index of the protare in the URR_protareInfos list. If negative, not in list. */
        bool m_hasURR_probabilityTables;                                            /**< *true* if URR probability tables present and *false* otherwise. */
        double m_URR_domainMin;                                                     /**< If URR probability tables present this is the minimum of the projectile energy domain for the tables. */
        double m_URR_domainMax;                                                     /**< If URR probability tables present this is the maximum of the projectile energy domain for the tables. */
        Vector<double> m_projectileMultiGroupBoundaries;                            /**< The multi-group boundaries for the projectile. Only used if m_crossSectionLookupMode and/or m_other1dDataLookupMode is multiGroup. */
        Vector<double> m_projectileMultiGroupBoundariesCollapsed;                   /**< The collased, multi-group boundaries for the projectile. Only used if m_crossSectionLookupMode and/or m_other1dDataLookupMode is multiGroup. */ 
        Vector<double> m_upscatterModelAGroupVelocities;                            /**< The speed of the projectile at each multi-group boundary. Need by upscatter model A. */

        Vector<Reaction *> m_reactions;                                             /**< The list of reactions. */
        Vector<Reaction *> m_orphanProducts;                                        /**< The list of orphan products. */
        bool m_isPhotoAtomic;                                                       /**< *true* if photo-atomic protare and false otherwise. */
        bool m_continuousEnergy;                                                    /**< If *true*, protare has continuous energy cross sections; otherwise, multi-group cross sections. */
        bool m_fixedGrid;                                                           /**< If *true*, continuous energy cross sections are fixed grid. */
        HeatedCrossSectionsContinuousEnergy m_heatedCrossSections;                  /**< Stores all cross section data for total and all reactions for all requested temperatures. */
        HeatedCrossSectionsMultiGroup m_heatedMultigroupCrossSections;              /**< Stores all multi-group cross section data for total and all reactions for all requested temperatures. */

        Vector<NuclideGammaBranchStateInfo *> m_nuclideGammaBranchStateInfos;       /**< List of all gamma branches for a nuclide. */
        Vector<NuclideGammaBranchInfo *> m_branches;                                /**< Condensed data on a nuclide's gamma branch including the gamma's energy, probability and the nuclide's residual state. */

        LUPI_HOST void setupNuclideGammaBranchStateInfos( SetupInfo &a_setupInfo, GIDI::ProtareSingle const &a_protare,
                bool a_makePhotonEmissionProbabilitiesOne, bool a_zeroNuclearLevelEnergyWidth );

    public:
        LUPI_HOST_DEVICE ProtareSingle( );
        LUPI_HOST ProtareSingle( LUPI::StatusMessageReporting &a_smr, GIDI::ProtareSingle const &a_protare, PoPI::Database const &a_pops, Transporting::MC &a_settings, 
                GIDI::Transporting::Particles const &a_particles, DomainHash const &a_domainHash, GIDI::Styles::TemperatureInfos const &a_temperatureInfos,
                std::set<int> const &a_reactionsToExclude, int a_reactionsToExcludeOffset = 0, bool a_allowFixedGrid = true );
        LUPI_HOST_DEVICE ~ProtareSingle( );

        LUPI_HOST_DEVICE bool isPhotoAtomic( ) const { return( m_isPhotoAtomic ); }
        LUPI_HOST_DEVICE bool continuousEnergy( ) const { return( m_continuousEnergy ); }
        LUPI_HOST_DEVICE bool fixedGrid( ) const { return( m_fixedGrid ); }
        LUPI_HOST_DEVICE HeatedCrossSectionsContinuousEnergy const &heatedCrossSections( ) const { return( m_heatedCrossSections ); }  /**< Returns a reference to the **m_heatedCrossSections** member. */
        LUPI_HOST_DEVICE HeatedCrossSectionsContinuousEnergy &heatedCrossSections( ) { return( m_heatedCrossSections ); }              /**< Returns a reference to the **m_heatedCrossSections** member. */
        LUPI_HOST_DEVICE HeatedCrossSectionsMultiGroup const &heatedMultigroupCrossSections( ) const { return( m_heatedMultigroupCrossSections ); } /**< Returns a reference to the **m_heatedMultigroupCrossSections** member. */
        LUPI_HOST_DEVICE HeatedCrossSectionsMultiGroup &heatedMultigroupCrossSections( ) { return( m_heatedMultigroupCrossSections ); } /**< Returns a reference to the **m_heatedMultigroupCrossSections** member. */

        LUPI_HOST_DEVICE const Vector<NuclideGammaBranchStateInfo *> &nuclideGammaBranchStateInfos( ) const { return( m_nuclideGammaBranchStateInfos ); }
                                                                                    /**< Returns a reference to the **m_nuclideGammaBranchStateInfos** member. */
        LUPI_HOST_DEVICE const Vector<NuclideGammaBranchInfo *> &branches( ) const { return( m_branches ); }
                                                                                    /**< Returns a reference to the **m_branches** member. */

// FIXME (1) see FIXME (1) in MC class.
        LUPI_HOST_DEVICE Vector<Reaction *> const &reactions( ) const { return( m_reactions ); }                 /**< Returns the value of the **m_reactions** member. */
// FIXME (1) see FIXME (1) in MC class.
        LUPI_HOST_DEVICE Vector<Reaction *> const &orphanProducts( ) const { return( m_orphanProducts ); }       /**< Returns the value of the **m_orphanProducts** member. */

        template <typename RNG, typename PUSHBACK>
        inline LUPI_HOST_DEVICE void sampleBranchingGammas( Sampling::Input &a_input, double a_projectileEnergy, int a_initialStateIndex, 
                RNG && a_rng, PUSHBACK && push_back, Sampling::ProductHandler &a_products ) const ;
        LUPI_HOST void setUserParticleIndex2( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid2( int a_particleIntid, int a_userParticleIndex );

// The rest are virtual methods defined in the Protare class.

        LUPI_HOST_DEVICE std::size_t numberOfProtares( ) const { return( 1 ); }                        /**< Returns the number of protares contained in *this*. */
        LUPI_HOST_DEVICE ProtareSingle const *protare( std::size_t a_index ) const ;
        LUPI_HOST_DEVICE ProtareSingle       *protare( std::size_t a_index );
        LUPI_HOST_DEVICE ProtareSingle const *protareWithReaction( int a_index ) const ;

        LUPI_HOST_DEVICE double minimumEnergy( ) const { 
            if( m_continuousEnergy ) return( m_heatedCrossSections.minimumEnergy( ) );
            return( m_heatedMultigroupCrossSections.minimumEnergy( ) ); }                                   /**< Returns the minimum cross section domain. */
        LUPI_HOST_DEVICE double maximumEnergy( ) const { 
            if( m_continuousEnergy ) return( m_heatedCrossSections.maximumEnergy( ) );
            return( m_heatedMultigroupCrossSections.maximumEnergy( ) ); }                                   /**< Returns the maximum cross section domain. */
        LUPI_HOST_DEVICE Vector<double> temperatures( std::size_t a_index = 0 ) const ;

        LUPI_HOST Vector<double> const &projectileMultiGroupBoundaries( ) const { return( m_projectileMultiGroupBoundaries ); }
                                                                                                            /**< Returns the value of the **m_projectileMultiGroupBoundaries** member. */
        LUPI_HOST Vector<double> const &projectileMultiGroupBoundariesCollapsed( ) const { return( m_projectileMultiGroupBoundariesCollapsed ); }
                                                                                                            /**< Returns the value of the **m_projectileMultiGroupBoundariesCollapsed** member. */

        LUPI_HOST_DEVICE std::size_t numberOfReactions( ) const { return( m_reactions.size( ) ); }                       /**< Returns the number of reactions of *this*. */
        LUPI_HOST_DEVICE Reaction const *reaction( int a_index ) const { return( m_reactions[a_index] ); }               /**< Returns the (a_index-1)^th reaction of *this*. */
        LUPI_HOST_DEVICE std::size_t numberOfOrphanProducts( ) const { return( m_orphanProducts.size( ) ); }             /**< Returns the number of orphan products of *this*. */
        LUPI_HOST_DEVICE Reaction const *orphanProduct( int a_index ) const { return( m_orphanProducts[a_index] ); }     /**< Returns the (a_index-1)^th orphan product of *this*. */

        LUPI_HOST_DEVICE bool hasFission( ) const ;
        LUPI_HOST_DEVICE String interaction( ) const { return( m_interaction ); }
        LUPI_HOST_DEVICE bool hasIncoherentDoppler( ) const ;

        LUPI_HOST_DEVICE int URR_index( ) const { return( m_URR_index ); }
        LUPI_HOST_DEVICE void URR_index( int a_URR_index ) { m_URR_index = a_URR_index; }
        LUPI_HOST_DEVICE bool inURR( double a_energy ) const ;
        LUPI_HOST_DEVICE bool hasURR_probabilityTables( ) const { return( m_hasURR_probabilityTables ); }
        LUPI_HOST_DEVICE double URR_domainMin( ) const { return( m_URR_domainMin ); }
        LUPI_HOST_DEVICE double URR_domainMax( ) const { return( m_URR_domainMax ); }
        LUPI_HOST_DEVICE bool reactionHasURR_probabilityTables( int a_index ) const { return( m_heatedCrossSections.reactionHasURR_probabilityTables( a_index ) ); }

        LUPI_HOST_DEVICE double threshold( std::size_t a_index ) const {
            if( m_continuousEnergy ) return( m_heatedCrossSections.threshold( a_index ) );
            return( m_heatedMultigroupCrossSections.threshold( a_index ) ); }                                       /**< Returns the threshold for the reaction at index *a_index*. */

        LUPI_HOST_DEVICE double crossSection(                              URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE void crossSectionVector( double a_temperature, double a_userFactor, int a_numberAllocated, double *a_crossSectionVector ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos,                  double a_temperature, double a_energy ) const ;
        template <typename RNG>
        inline LUPI_HOST_DEVICE int sampleReaction(                               URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, double a_crossSection, RNG && a_rng  ) const ;

        LUPI_HOST_DEVICE double depositionEnergy(   int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double depositionMomentum( int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double productionEnergy(   int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double gain(               int a_hashIndex, double a_temperature, double a_energy, int a_particleIndex ) const ;
        LUPI_HOST_DEVICE double gainViaIntid(       int a_hashIndex, double a_temperature, double a_energy, int a_particleIntid ) const ;

        LUPI_HOST_DEVICE Vector<double> const &upscatterModelAGroupVelocities( ) const { return( m_upscatterModelAGroupVelocities ); }   /**< Returns a reference to the **m_upscatterModelAGroupVelocities** member. */

        LUPI_HOST_DEVICE void serialize2( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST_DEVICE long sizeOf2( ) const { return sizeof(*this); }
};

/*
============================================================
===================== ProtareComposite =====================
============================================================
*/
class ProtareComposite : public Protare {

    private:
        Vector<ProtareSingle *> m_protares;                              /**< List of protares added to *this* instance. */
        std::size_t m_numberOfReactions;                                    /**< The sum of the number of reaction for all stored protares. */
        std::size_t m_numberOfOrphanProducts;                               /**< The sum of the number of reaction for all stored protares. */
        double m_minimumEnergy;                                             /**< The maximum of the minimum cross section domains. */
        double m_maximumEnergy;                                             /**< The minimum of the maximum cross section domains. */

    public:
        LUPI_HOST_DEVICE ProtareComposite( );
        LUPI_HOST ProtareComposite( LUPI::StatusMessageReporting &a_smr, GIDI::ProtareComposite const &a_protare, PoPI::Database const &a_pops, Transporting::MC &a_settings, 
                GIDI::Transporting::Particles const &a_particles, DomainHash const &a_domainHash, GIDI::Styles::TemperatureInfos const &a_temperatureInfos,
                std::set<int> const &a_reactionsToExclude, int a_reactionsToExcludeOffset = 0, bool a_allowFixedGrid = true );
        LUPI_HOST_DEVICE ~ProtareComposite( );

        Vector<ProtareSingle *> protares( ) const { return( m_protares ); }       /**< Returns the value of the **m_protares** member. */
        LUPI_HOST void setUserParticleIndex2( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid2( int a_particleIntid, int a_userParticleIndex );

// The rest are virtual methods defined in the Protare class.

        LUPI_HOST_DEVICE std::size_t numberOfProtares( ) const { return( m_protares.size( ) ); }     /**< Returns the number of protares contained in *this*. */
        LUPI_HOST_DEVICE ProtareSingle const *protare( std::size_t a_index ) const ;
        LUPI_HOST_DEVICE ProtareSingle       *protare( std::size_t a_index );
        LUPI_HOST_DEVICE ProtareSingle const *protareWithReaction( int a_index ) const ;

        LUPI_HOST_DEVICE double minimumEnergy( ) const { return( m_minimumEnergy ); }     /**< Returns the value of the **m_minimumEnergy** member. */
        LUPI_HOST_DEVICE double maximumEnergy( ) const { return( m_maximumEnergy ); }     /**< Returns the value of the **m_maximumEnergy** member. */
        LUPI_HOST_DEVICE Vector<double> temperatures( std::size_t a_index = 0 ) const ;

        LUPI_HOST Vector<double> const &projectileMultiGroupBoundaries( ) const { return( m_protares[0]->projectileMultiGroupBoundaries( ) ); }    
                                                                            /**< Returns the value of the **m_projectileMultiGroupBoundaries** member. */
        LUPI_HOST Vector<double> const &projectileMultiGroupBoundariesCollapsed( ) const { return( m_protares[0]->projectileMultiGroupBoundariesCollapsed( ) ); }
                                                                            /**< Returns the value of the **m_projectileMultiGroupBoundariesCollapsed** member. */

        LUPI_HOST_DEVICE std::size_t numberOfReactions( ) const { return( m_numberOfReactions ); }
                                                                            /**< Returns the value of the **m_numberOfReactions** member. */
        LUPI_HOST_DEVICE Reaction const *reaction( int a_index ) const ;
        LUPI_HOST_DEVICE std::size_t numberOfOrphanProducts( ) const { return( m_numberOfOrphanProducts ); }
                                                                            /**< Returns the value of the **m_numberOfOrphanProducts** member. */
        LUPI_HOST_DEVICE Reaction const *orphanProduct( int a_index ) const ;

        LUPI_HOST_DEVICE bool hasFission( ) const ;
        LUPI_HOST_DEVICE bool hasIncoherentDoppler( ) const ;

        LUPI_HOST_DEVICE int URR_index( ) const { return( -1 ); }
        LUPI_HOST_DEVICE bool hasURR_probabilityTables( ) const ;
        LUPI_HOST_DEVICE double URR_domainMin( ) const ;
        LUPI_HOST_DEVICE double URR_domainMax( ) const ;
        LUPI_HOST_DEVICE bool reactionHasURR_probabilityTables( int a_index ) const ;

        LUPI_HOST_DEVICE double threshold( std::size_t a_index ) const ;

        LUPI_HOST_DEVICE double crossSection(                              URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE void crossSectionVector( double a_temperature, double a_userFactor, int a_numberAllocated, double *a_crossSectionVector ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos,                  double a_temperature, double a_energy ) const ;
        template <typename RNG>
        inline LUPI_HOST_DEVICE int sampleReaction(                               URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, double a_crossSection, RNG && a_rng ) const ;

        LUPI_HOST_DEVICE double depositionEnergy(   int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double depositionMomentum( int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double productionEnergy(   int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double gain(               int a_hashIndex, double a_temperature, double a_energy, int a_particleIndex ) const ;
        LUPI_HOST_DEVICE double gainViaIntid(       int a_hashIndex, double a_temperature, double a_energy, int a_particleIntid ) const ;

        LUPI_HOST_DEVICE Vector<double> const &upscatterModelAGroupVelocities( ) const { return( m_protares[0]->upscatterModelAGroupVelocities( ) ); }
                                                                            /**< Returns a reference to the **m_upscatterModelAGroupVelocities** member. */

        LUPI_HOST_DEVICE void serialize2( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST_DEVICE long sizeOf2( ) const { return sizeof(*this); }
};

/*
============================================================
======================== ProtareTNSL =======================
============================================================
*/
class ProtareTNSL : public Protare {

    private:
        std::size_t m_numberOfTNSLReactions;                                /**< The number of reactions of the TNSL protare. */
        double m_TNSL_maximumEnergy;                                        /**< The maximum energy of the cross section domain for the TNSL protare. */
        double m_TNSL_maximumTemperature;                                   /**< The highest temperature for processed data for the TNSL protare. */
        ProtareSingle *m_protareWithElastic;                                /**< Protare with non thermal neutron scattering law data. */
        ProtareSingle *m_TNSL;                                              /**< Protare with thermal neutron scattering law data. */
        ProtareSingle *m_protareWithoutElastic;                             /**< Same as *m_protare* but without elastic. */

    public:
        LUPI_HOST_DEVICE ProtareTNSL( );
        LUPI_HOST ProtareTNSL( LUPI::StatusMessageReporting &a_smr, GIDI::ProtareTNSL const &a_protare, PoPI::Database const &a_pops, Transporting::MC &a_settings, 
                GIDI::Transporting::Particles const &a_particles, DomainHash const &a_domainHash, GIDI::Styles::TemperatureInfos const &a_temperatureInfos,
                std::set<int> const &a_reactionsToExclude, int a_reactionsToExcludeOffset = 0, bool a_allowFixedGrid = true );
        LUPI_HOST_DEVICE ~ProtareTNSL( );

        LUPI_HOST_DEVICE ProtareSingle const *protareWithElastic( ) const { return( m_protareWithElastic ); }        /**< Returns the **m_protareWithElastic** member. */
        LUPI_HOST_DEVICE ProtareSingle const *TNSL( ) const { return( m_TNSL ); }                                    /**< Returns the **m_TNSL** member. */
        LUPI_HOST_DEVICE ProtareSingle const *protareWithoutElastic( ) const { return( m_protareWithoutElastic ); }  /**< Returns the **m_protareWithoutElastic** member. */

        LUPI_HOST_DEVICE double TNSL_maximumEnergy( ) const { return( m_TNSL_maximumEnergy ); }
        LUPI_HOST_DEVICE double TNSL_maximumTemperature( ) const { return( m_TNSL_maximumTemperature ); }
        LUPI_HOST void setUserParticleIndex2( int a_particleIndex, int a_userParticleIndex );
        LUPI_HOST void setUserParticleIndexViaIntid2( int a_particleIntid, int a_userParticleIndex );

// The rest are virtual methods defined in the Protare class.

        LUPI_HOST_DEVICE std::size_t numberOfProtares( ) const { return( 2 ); }  /**< Always Returns 2. */
        LUPI_HOST_DEVICE ProtareSingle const *protare( std::size_t a_index ) const ;
        LUPI_HOST_DEVICE ProtareSingle       *protare( std::size_t a_index );
        LUPI_HOST_DEVICE ProtareSingle const *protareWithReaction( int a_index ) const ;

        LUPI_HOST_DEVICE double minimumEnergy( ) const { return( m_protareWithElastic->minimumEnergy( ) ); }   /**< Returns the minimum cross section domain. */
        LUPI_HOST_DEVICE double maximumEnergy( ) const { return( m_protareWithElastic->maximumEnergy( ) ); }   /**< Returns the maximum cross section domain. */
        LUPI_HOST_DEVICE Vector<double> temperatures( std::size_t a_index = 0 ) const ;

        LUPI_HOST Vector<double> const &projectileMultiGroupBoundaries( ) const { return( m_protareWithElastic->projectileMultiGroupBoundaries( ) ); }
                                                                            /**< Returns the value of the **m_projectileMultiGroupBoundaries** member. */
        LUPI_HOST Vector<double> const &projectileMultiGroupBoundariesCollapsed( ) const { return( m_protareWithElastic->projectileMultiGroupBoundariesCollapsed( ) ); }
                                                                            /**< Returns the value of the **m_projectileMultiGroupBoundariesCollapsed** member. */

        LUPI_HOST_DEVICE std::size_t numberOfReactions( ) const { return( m_TNSL->numberOfReactions( ) + m_protareWithElastic->numberOfReactions( ) ); }
        LUPI_HOST_DEVICE Reaction const *reaction( int a_index ) const ;
        LUPI_HOST_DEVICE std::size_t numberOfOrphanProducts( ) const { return( m_protareWithElastic->numberOfOrphanProducts( ) ); }
                                                                            /**< Returns the number of orphan products in the normal ProtareSingle. */
        LUPI_HOST_DEVICE Reaction const *orphanProduct( int a_index ) const { return( m_protareWithElastic->orphanProduct( a_index ) ); }
                                                                            /**< Returns the (a_index - 1 )^th orphan product in the normal ProtareSingle. */

        LUPI_HOST_DEVICE bool hasFission( ) const { return( m_protareWithElastic->hasFission( ) ); }    /* Returns the normal ProtareSingle's hasFission value. */
        LUPI_HOST_DEVICE bool hasIncoherentDoppler( ) const { return( false ); }                        /* Always returns false as this is a neutron as projectile and not a photon. */

        LUPI_HOST_DEVICE int URR_index( ) const { return( -1 ); }
        LUPI_HOST_DEVICE bool hasURR_probabilityTables( ) const { return( m_protareWithElastic->hasURR_probabilityTables( ) ); }
        LUPI_HOST_DEVICE double URR_domainMin( ) const { return( m_protareWithElastic->URR_domainMin( ) ); }
        LUPI_HOST_DEVICE double URR_domainMax( ) const { return( m_protareWithElastic->URR_domainMax( ) ); }
        LUPI_HOST_DEVICE bool reactionHasURR_probabilityTables( int a_index ) const ;

        LUPI_HOST_DEVICE double threshold( std::size_t a_index ) const ;

        LUPI_HOST_DEVICE double crossSection(                              URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE void crossSectionVector( double a_temperature, double a_userFactor, int a_numberAllocated, double *a_crossSectionVector ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, bool a_sampling = false ) const ;
        LUPI_HOST_DEVICE double reactionCrossSection( int a_reactionIndex, URR_protareInfos const &a_URR_protareInfos,                  double a_temperature, double a_energy ) const ;
        template <typename RNG>
        inline LUPI_HOST_DEVICE int sampleReaction(                               URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, double a_crossSection, RNG && a_rng ) const ;

        LUPI_HOST_DEVICE double depositionEnergy(   int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double depositionMomentum( int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double productionEnergy(   int a_hashIndex, double a_temperature, double a_energy ) const ;
        LUPI_HOST_DEVICE double gain(               int a_hashIndex, double a_temperature, double a_energy, int a_particleIndex ) const ;
        LUPI_HOST_DEVICE double gainViaIntid(       int a_hashIndex, double a_temperature, double a_energy, int a_particleIntid ) const ;

        LUPI_HOST_DEVICE Vector<double> const &upscatterModelAGroupVelocities( ) const { return( m_protareWithElastic->upscatterModelAGroupVelocities( ) ); }
                                                                            /**< Returns a reference to the **m_upscatterModelAGroupVelocities** member. */

        LUPI_HOST_DEVICE void serialize2( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST_DEVICE long sizeOf2( ) const { return sizeof(*this); }
};

/*
============================================================
=========================== Others =========================
============================================================
*/
LUPI_HOST Protare *protareFromGIDIProtare( LUPI::StatusMessageReporting &a_smr, GIDI::Protare const &a_protare, PoPI::Database const &a_pops, Transporting::MC &a_settings, GIDI::Transporting::Particles const &a_particles,
                DomainHash const &a_domainHash, GIDI::Styles::TemperatureInfos const &a_temperatureInfos, std::set<int> const &a_reactionsToExclude,
                int a_reactionsToExcludeOffset = 0, bool a_allowFixedGrid = true );
LUPI_HOST Vector<double> GIDI_VectorDoublesToMCGIDI_VectorDoubles( GIDI::Vector a_vector );
LUPI_HOST void addVectorItemsToSet( Vector<int> const &a_from, std::set<int> &a_to );

LUPI_HOST_DEVICE int distributionTypeToInt( Distributions::Type a_type );
LUPI_HOST_DEVICE Distributions::Type intToDistributionType( int a_type );
LUPI_HOST_DEVICE void serializeProducts( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Vector<Product *> &a_products );
LUPI_HOST_DEVICE void serializeDelayedNeutrons( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Vector<DelayedNeutron *> &a_delayedNeutrons );
LUPI_HOST_DEVICE void serializeQs( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Vector<Functions::Function1d_d1 *> &a_Qs );
LUPI_HOST_DEVICE void serializeFissionResiduals( GIDI::Construction::FissionResiduals &a_fissionResiduals, 
                LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );

LUPI_HOST void convertACE_URR_probabilityTablesFromGIDI( GIDI::ProtareSingle const &a_protare, Transporting::MC &a_settings, SetupInfo &a_setupInfo );
LUPI_HOST_DEVICE Transporting::URR_mode serializeURR_mode( Transporting::URR_mode a_URR_mode, LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
LUPI_HOST_DEVICE ACE_URR_probabilityTables *serializeACE_URR_probabilityTables( ACE_URR_probabilityTables *a_ACE_URR_probabilityTables,
                LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );

LUPI_HOST_DEVICE Distributions::Distribution *serializeDistribution( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, 
                Distributions::Distribution *a_distribution );

LUPI_HOST std::vector<double> vectorToSTD_vector( Vector<double> a_input );
LUPI_HOST std::vector<double> vectorToSTD_vector( Vector<float> a_input );

}           // End of namespace MCGIDI.

#include "MCGIDI_headerSource.hpp"

#endif      // End of MCGIDI_hpp_included
