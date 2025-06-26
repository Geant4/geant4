/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

namespace DoubleDifferentialCrossSection {

/* *********************************************************************************************************//**
 * Returns the name of the DebyeWallerIntegral child in *a_node*, independent of the GNDS format.
 *
 * @param a_node            [in]    The **HAPI::Node** whose child are checked.
 ***********************************************************************************************************/

static char const *getDebyeWallerIntegralName( HAPI::Node const &a_node ) {

    HAPI::Node const &node = a_node.child( GIDI_DebyeWallerIntegralChars );
    if( node.empty( ) ) return( GIDI_DebyeWallerChars );

    return( GIDI_DebyeWallerIntegralChars );
}

/* *********************************************************************************************************//**
 * Returns the boundAtomCrossSection child in *a_node*, independent of the GNDS format.
 *
 * @param a_node            [in]    The **HAPI::Node** whose child are checked.
 ***********************************************************************************************************/

static char const *getBoundAtomCrossSectionName( HAPI::Node const &a_node ) {

    HAPI::Node const &node = a_node.child( GIDI_boundAtomCrossSectionChars );
    if( node.empty( ) ) return( GIDI_characteristicCrossSectionChars );

    return( GIDI_boundAtomCrossSectionChars );
}

/*! \class Base
 * Base class inherited by DoubleDifferentialCrossSection forms.
 */

/* *********************************************************************************************************//**
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the Base.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_type            [in]    The FormType for the DoubleDifferentialCrossSection form.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Base::Base( HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type, Suite *a_parent ) :
        Form( a_node, a_setupInfo, a_type, a_parent ) {

}

/*! \class CoherentPhotoAtomicScattering
 * This is the **coherentPhotonScattering** style class.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering( Construction::Settings const &a_construction, HAPI::Node const &a_node,
		SetupInfo &a_setupInfo, LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, Suite *a_parent ) :
        Base( a_node, a_setupInfo, FormType::coherentPhotonScattering, a_parent ),
        m_formFactor( data1dParse( a_construction, a_node.child( GIDI_formFactorChars ).first_child( ), a_setupInfo, nullptr ) ),
        m_realAnomalousFactor( data1dParseAllowEmpty( a_construction, a_node.child( GIDI_realAnomalousFactorChars ).first_child( ), a_setupInfo, nullptr ) ),
        m_imaginaryAnomalousFactor( data1dParseAllowEmpty( a_construction, a_node.child( GIDI_imaginaryAnomalousFactorChars ).first_child( ), a_setupInfo, nullptr ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

CoherentPhotoAtomicScattering::~CoherentPhotoAtomicScattering( ) {

    delete m_formFactor;
    delete m_realAnomalousFactor;
    delete m_imaginaryAnomalousFactor;
}

/*! \class IncoherentPhotoAtomicScattering
 * This is the **incoherentPhotonScattering** class.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

IncoherentPhotoAtomicScattering::IncoherentPhotoAtomicScattering( Construction::Settings const &a_construction, HAPI::Node const &a_node,
		SetupInfo &a_setupInfo, LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, Suite *a_parent ) :
        Base( a_node, a_setupInfo, FormType::incoherentPhotonScattering, a_parent ),
        m_scatteringFactor( nullptr ) {

    HAPI::Node const scatteringFactorChild = a_node.child( GIDI_scatteringFactorChars );
    if( scatteringFactorChild.empty( ) ) {
        m_scatteringFactor = data1dParse( a_construction, a_node.first_child( ), a_setupInfo, nullptr ); }
    else {
        m_scatteringFactor = data1dParse( a_construction, scatteringFactorChild.first_child( ), a_setupInfo, nullptr );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

IncoherentPhotoAtomicScattering::~IncoherentPhotoAtomicScattering( ) {

    delete m_scatteringFactor;
}

/*! \class IncoherentPhotoAtomicScattering
 * This is the **incoherentPhotonScattering** class.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

IncoherentBoundToFreePhotoAtomicScattering::IncoherentBoundToFreePhotoAtomicScattering( Construction::Settings const &a_construction, HAPI::Node const &a_node,
		SetupInfo &a_setupInfo, LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, Suite *a_parent ) :
        Base( a_node, a_setupInfo, FormType::incoherentBoundToFreePhotonScattering, a_parent ),
        m_ComptonProfile( nullptr ) {

    HAPI::Node const ComptonProfileChild = a_node.child( GIDI_ComptonProfileChars );
    m_ComptonProfile = data1dParse( a_construction, ComptonProfileChild.first_child( ), a_setupInfo, nullptr );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

IncoherentBoundToFreePhotoAtomicScattering::~IncoherentBoundToFreePhotoAtomicScattering( ) {

    delete m_ComptonProfile;
}

namespace n_ThermalNeutronScatteringLaw {

/*! \class S_table
 * This class represents a **GNDS** cumulative scattering factor **S_table** instance which is a function of temperaure 
 * \f$T\f$ and energy \f$E\f$ as \f$S(T,E)\f$.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

S_table::S_table( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::generic, nullptr ),
        m_function2d( data2dParse( a_construction, a_node.first_child( ), a_setupInfo, nullptr ) ) {

    m_function2d->setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

S_table::~S_table( ) {

    delete m_function2d;

}

/*! \class CoherentElastic
 * This is the **CoherentElastic** class.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the **PoPs** node under the **reactionSuite** node.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

CoherentElastic::CoherentElastic( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, Suite *a_parent ) :
        Base( a_node, a_setupInfo, FormType::coherentElastic, a_parent ),
        m_S_table( a_construction, a_node.child( GIDI_S_tableChars ), a_setupInfo ) {

    m_S_table.setAncestor( this );        
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

CoherentElastic::~CoherentElastic( ) {

}

/*! \class DebyeWallerIntegral
 * This class represents a **GNDS** Debye-Waller integral **DebyeWallerIntegral** which is a function of temperaure \f$T\f$ as \f$W'(T)\f$.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

DebyeWallerIntegral::DebyeWallerIntegral( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::generic, nullptr ),
        m_function1d( data1dParse( a_construction, a_node.first_child( ), a_setupInfo, nullptr ) ) {

    m_function1d->setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

DebyeWallerIntegral::~DebyeWallerIntegral( ) {

    delete m_function1d;
}

/*! \class IncoherentElastic
 * This is the **IncoherentElastic** class.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

IncoherentElastic::IncoherentElastic( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		LUPI_maybeUnused PoPI::Database const &a_pops, LUPI_maybeUnused PoPI::Database const &a_internalPoPs, Suite *a_parent ) :
        Base( a_node, a_setupInfo, FormType::incoherentElastic, a_parent ),
        m_boundAtomCrossSection( a_node.child( getBoundAtomCrossSectionName( a_node ) ), a_setupInfo ),
        m_DebyeWallerIntegral( a_construction, a_node.child( getDebyeWallerIntegralName( a_node ) ), a_setupInfo ) {

    m_boundAtomCrossSection.setAncestor( this );        
    m_DebyeWallerIntegral.setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

IncoherentElastic::~IncoherentElastic( ) {

}

/*! \class Options
 * This is the **options** class.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

Options::Options( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::generic, nullptr ),
        m_calculatedAtThermal( strcmp( a_node.attribute_as_string( GIDI_calculatedAtThermalChars ).c_str( ), GIDI_trueChars ) == 0 ),
        m_asymmetric( strcmp( a_node.attribute_as_string( GIDI_asymmetricChars ).c_str( ), GIDI_trueChars ) == 0 ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Options::~Options( ) {

}

/*! \class T_effective
 * This class represents a **GNDS** **T_effective** which is a function of temperaure \f$T_{\rm eff}(T)\f$.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

T_effective::T_effective( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::generic, nullptr ),
        m_function1d( data1dParseAllowEmpty( a_construction, a_node.first_child( ), a_setupInfo, nullptr ) ) {

    if( m_function1d != nullptr ) m_function1d->setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

T_effective::~T_effective( ) {

    delete m_function1d;
}

/*! \class ScatteringAtom
 * This is the **scatteringAtom** class.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

ScatteringAtom::ScatteringAtom( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::generic, nullptr ),
        m_mass( a_node.child( GIDI_massChars ), a_setupInfo ),
        m_freeAtomCrossSection( a_node.child( GIDI_freeAtomCrossSectionChars ), a_setupInfo ), 
        m_e_critical( a_node.child( GIDI_e_criticalChars ), a_setupInfo ),
        m_e_max( a_node.child( GIDI_e_maxChars ), a_setupInfo ),
        m_T_effective( a_construction, a_node.child( GIDI_T_effectiveChars ), a_setupInfo ) {

    m_mass.setAncestor( this );
    m_freeAtomCrossSection.setAncestor( this );
    m_e_critical.setAncestor( this );
    m_e_max.setAncestor( this );
    m_T_effective.setAncestor( this );

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

ScatteringAtom::~ScatteringAtom( ) {

}

/*! \class S_alpha_beta
 * This class represents a **GNDS** **S_alpha_beta** which is \f$S(T,\alpha,\beta)\f$.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 ***********************************************************************************************************/

S_alpha_beta::S_alpha_beta( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo ) :
        Form( a_node, a_setupInfo, FormType::generic, nullptr ),
        m_function3d( nullptr ) { // data3dParse( a_construction, a_node.first_child( ), a_setupInfo, nullptr ) )

// FIXME BRB
//    m_function3d->setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

S_alpha_beta::~S_alpha_beta( ) {

    delete m_function3d;
}

/*! \class IncoherentInelastic
 * This is the **IncoherentInelastic** class.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options for parsing.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

IncoherentInelastic::IncoherentInelastic( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Suite *a_parent ) :
        Base( a_node, a_setupInfo, FormType::incoherentInelastic, a_parent ),
        m_options( a_construction, a_node.child( GIDI_optionsChars ), a_setupInfo ),
        m_scatteringAtoms( a_construction, GIDI_scatteringAtomsChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseScatteringAtom, nullptr ),
        m_S_alpha_beta( a_construction, a_node.child( GIDI_S_alpha_betaChars ), a_setupInfo ) {

    m_options.setAncestor( this );
    m_scatteringAtoms.setAncestor( this );
    m_S_alpha_beta.setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

IncoherentInelastic::~IncoherentInelastic( ) {

}

}                   // End namespace n_ThermalNeutronScatteringLaw.

}                   // End namespace DoubleDifferentialCrossSection.

}                   // End namespace GIDI.
