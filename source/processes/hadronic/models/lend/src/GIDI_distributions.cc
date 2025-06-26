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

namespace Distributions {

/*! \class Distribution
 * This is the base class inherited by all distribution form classes.
 */

/* *********************************************************************************************************//**
 *
 * @param a_moniker         [in]    The **GNDS** moniker for the distribution.
 * @param a_type            [in]    The FormType for the distribution form.
 * @param a_label           [in]    The label for the distribution.
 * @param a_productFrame    [in]    The frame the product data are in.
 ***********************************************************************************************************/

Distribution::Distribution( std::string const &a_moniker, FormType a_type, std::string const &a_label, Frame a_productFrame ) :
        Form( a_moniker, a_type, a_label ),
        m_productFrame( a_productFrame ) {

}

/* *********************************************************************************************************//**
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the distribution.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_type            [in]    The FormType for the distribution form.
 * @param a_parent          [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

Distribution::Distribution( HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type, Suite *a_parent ) :
        Form( a_node, a_setupInfo, a_type, a_parent ),
        m_productFrame( parseFrame( a_node, a_setupInfo, GIDI_productFrameChars ) ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Distribution::toXMLNodeStarter( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;

    attributes += a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
    attributes += a_writeInfo.addAttribute( GIDI_productFrameChars, frameToString( m_productFrame ) );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
}

/*! \class MultiGroup3d
 * Class for the GNDS **multiGroup3d** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the MultiGroup3d.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

MultiGroup3d::MultiGroup3d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::multiGroup3d, a_parent ),
        m_gridded3d( a_construction, a_node.child( GIDI_gridded3dChars ), a_setupInfo ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void MultiGroup3d::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    toXMLNodeStarter( a_writeInfo, a_indent );
    m_gridded3d.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class AngularTwoBody
 * Class for the GNDS **angularTwoBody** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_label           [in]    The label for *this* form.
 * @param a_productFrame    [in]    The frame the product data as in.
 * @param a_angular         [in]    The 2-d functional angular representation.
 ***********************************************************************************************************/

AngularTwoBody::AngularTwoBody( std::string const &a_label, Frame a_productFrame, Functions::Function2dForm *a_angular ) :
        Distribution( GIDI_angularTwoBodyChars, FormType::angularTwoBody, a_label, a_productFrame ),
        m_angular( a_angular ) {

    if( a_angular != nullptr ) a_angular->setAncestor( this );
}

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the AngularTwoBody.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

AngularTwoBody::AngularTwoBody( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::angularTwoBody, a_parent ),
        m_angular( data2dParse( a_construction, a_node.first_child( ), a_setupInfo, nullptr ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

AngularTwoBody::~AngularTwoBody( ) {

    delete m_angular;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/
 
void AngularTwoBody::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    toXMLNodeStarter( a_writeInfo, a_indent );
    if( m_angular != nullptr ) m_angular->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class KalbachMann
 * Class for the GNDS **KalbachMann** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the KalbachMann.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

KalbachMann::KalbachMann( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::KalbachMann, a_parent ),
        m_f( data2dParse( a_construction, a_node.child( GIDI_fChars ).first_child( ), a_setupInfo, nullptr ) ),
        m_r( data2dParse( a_construction, a_node.child( GIDI_rChars ).first_child( ), a_setupInfo, nullptr ) ),
        m_a( nullptr ) {

    HAPI::Node const &aNode = a_node.child( GIDI_aChars );
    if( !aNode.empty( ) ) m_a = data2dParse( a_construction, aNode.first_child( ), a_setupInfo, nullptr );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

KalbachMann::~KalbachMann( ) {

    delete m_f;
    delete m_r;
    delete m_a;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/
 
void KalbachMann::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    toXMLNodeStarter( a_writeInfo, a_indent );
    m_f->toXMLList( a_writeInfo, indent2 );
    m_r->toXMLList( a_writeInfo, indent2 );
    if( m_a != nullptr ) m_a->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class EnergyAngular
 * Class for the GNDS **energyAngular** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the EnergyAngular.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

EnergyAngular::EnergyAngular( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::energyAngular, a_parent ),
        m_energyAngular( data3dParse( a_construction, a_node.first_child( ), a_setupInfo, nullptr ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

EnergyAngular::~EnergyAngular( ) {

    delete m_energyAngular;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/
 
void EnergyAngular::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    toXMLNodeStarter( a_writeInfo, a_indent );
    m_energyAngular->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class EnergyAngularMC
 * Class for the GNDS **energyAngularMC** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the EnergyAngularMC.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

EnergyAngularMC::EnergyAngularMC( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::energyAngularMC, a_parent ),
        m_energy( data2dParse( a_construction, a_node.child( GIDI_energyChars ).first_child( ), a_setupInfo, nullptr ) ),
        m_energyAngular( data3dParse( a_construction, a_node.child( GIDI_energyAngularChars ).first_child( ), a_setupInfo, nullptr ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

EnergyAngularMC::~EnergyAngularMC( ) {

    delete m_energy;
    delete m_energyAngular;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/
 
void EnergyAngularMC::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string indent3 = a_writeInfo.incrementalIndent( indent2 );

    toXMLNodeStarter( a_writeInfo, a_indent );

    a_writeInfo.addNodeStarter( indent2, GIDI_energyChars );
    m_energy->toXMLList( a_writeInfo, indent3 );
    a_writeInfo.addNodeEnder( GIDI_energyChars );

    a_writeInfo.addNodeStarter( indent2, GIDI_energyAngularChars );
    m_energyAngular->toXMLList( a_writeInfo, indent3 );
    a_writeInfo.addNodeEnder( GIDI_energyAngularChars );

    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class AngularEnergy
 * Class for the GNDS **angularEnergy** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the AngularEnergy.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

AngularEnergy::AngularEnergy( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::angularEnergy, a_parent ),
        m_angularEnergy( data3dParse( a_construction, a_node.first_child( ), a_setupInfo, nullptr ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

AngularEnergy::~AngularEnergy( ) {

    delete m_angularEnergy;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/
 
void AngularEnergy::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    toXMLNodeStarter( a_writeInfo, a_indent );
    m_angularEnergy->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class AngularEnergyMC
 * Class for the GNDS **angularEnergyMC** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the AngularEnergyMC.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

AngularEnergyMC::AngularEnergyMC( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::angularEnergyMC, a_parent ),
        m_angular( data2dParse( a_construction, a_node.child( GIDI_angularChars ).first_child( ), a_setupInfo, nullptr ) ),
        m_angularEnergy( data3dParse( a_construction, a_node.child( GIDI_angularEnergyChars ).first_child( ), a_setupInfo, nullptr ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

AngularEnergyMC::~AngularEnergyMC( ) {

    delete m_angular;
    delete m_angularEnergy;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/
 
void AngularEnergyMC::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string indent3 = a_writeInfo.incrementalIndent( indent2 );

    toXMLNodeStarter( a_writeInfo, a_indent );
    a_writeInfo.addNodeStarter( indent2, GIDI_angularChars );
    m_angular->toXMLList( a_writeInfo, indent3 );
    a_writeInfo.addNodeEnder( GIDI_angularChars );
    a_writeInfo.addNodeStarter( indent2, GIDI_angularEnergyChars );
    m_angularEnergy->toXMLList( a_writeInfo, indent3 );
    a_writeInfo.addNodeEnder( GIDI_angularEnergyChars );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class Uncorrelated
 * Class for the GNDS **uncorrelated** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the Uncorrelated.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

Uncorrelated::Uncorrelated( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::uncorrelated, a_parent ),
        m_angular( data2dParse( a_construction, a_node.child( GIDI_angularChars ).first_child( ), a_setupInfo, nullptr ) ),
        m_energy( data2dParse( a_construction, a_node.child( GIDI_energyChars ).first_child( ), a_setupInfo, nullptr ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Uncorrelated::~Uncorrelated( ) {

    delete m_angular;
    delete m_energy;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/
 
void Uncorrelated::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string indent3 = a_writeInfo.incrementalIndent( indent2 );

    toXMLNodeStarter( a_writeInfo, a_indent );

    a_writeInfo.addNodeStarter( indent2, GIDI_angularChars, "" );
    m_angular->toXMLList( a_writeInfo, indent3 );
    a_writeInfo.addNodeEnder( GIDI_angularChars );

    a_writeInfo.addNodeStarter( indent2, GIDI_energyChars, "" );
    m_energy->toXMLList( a_writeInfo, indent3 );
    a_writeInfo.addNodeEnder( GIDI_energyChars );

    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class LLNLAngularEnergy
 * Class for the GNDS **LLNLAngularEnergy** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the LLNLAngularEnergy.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

LLNLAngularEnergy::LLNLAngularEnergy( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::LLNL_angularEnergy, a_parent ),
        m_angular( data2dParse( a_construction, a_node.child( GIDI_LLNLAngularOfAngularEnergyChars ).first_child( ), a_setupInfo, nullptr ) ),
        m_angularEnergy( data3dParse( a_construction, a_node.child( GIDI_LLNLAngularEnergyOfAngularEnergyChars ).first_child( ), a_setupInfo, nullptr ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LLNLAngularEnergy::~LLNLAngularEnergy( ) {

    delete m_angular;
    delete m_angularEnergy;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/
 
void LLNLAngularEnergy::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string indent3 = a_writeInfo.incrementalIndent( indent2 );

    toXMLNodeStarter( a_writeInfo, a_indent );
    a_writeInfo.addNodeStarter( indent2, GIDI_LLNLAngularOfAngularEnergyChars );
    m_angular->toXMLList( a_writeInfo, indent3 );
    a_writeInfo.addNodeEnder( GIDI_LLNLAngularOfAngularEnergyChars );
    a_writeInfo.addNodeStarter( indent2, GIDI_LLNLAngularEnergyOfAngularEnergyChars );
    m_angularEnergy->toXMLList( a_writeInfo, indent3 );
    a_writeInfo.addNodeEnder( GIDI_LLNLAngularEnergyOfAngularEnergyChars );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class CoherentPhotoAtomicScattering
 * Class for the GNDS **coherentPhotoAtomicScattering** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the CoherentPhotoAtomicScattering.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, 
                SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::coherentPhotonScattering, a_parent ),
        m_href( a_node.attribute_as_string( GIDI_hrefChars ) ) {

}

/*! \class IncoherentPhotoAtomicScattering
 * Class for the GNDS **incoherentPhotoAtomicScattering** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the IncoherentPhotoAtomicScattering.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

IncoherentPhotoAtomicScattering::IncoherentPhotoAtomicScattering( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node,
                SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::incoherentPhotonScattering, a_parent ),
        m_href( a_node.attribute_as_string( GIDI_hrefChars ) ) {

}

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the IncoherentPhotoAtomicScattering.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

IncoherentBoundToFreePhotoAtomicScattering::IncoherentBoundToFreePhotoAtomicScattering( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node,
                SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::incoherentBoundToFreePhotonScattering, a_parent ),
        m_href( a_node.attribute_as_string( GIDI_hrefChars ) ) {
}

/*! \class ThermalNeutronScatteringLaw
 * Class for the GNDS **thermalNeutronScatteringLaw** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the ThermalNeutronScatteringLaw.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

ThermalNeutronScatteringLaw::ThermalNeutronScatteringLaw( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node,
		SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::thermalNeutronScatteringLaw, a_parent ),
        m_href( a_node.attribute_as_string( GIDI_hrefChars ) ) {

}

/*! \class Branching3d
 * Class for the GNDS **branching3d** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the Branching3d.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

Branching3d::Branching3d( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::branching3d, a_parent ),
        m_initialState( a_setupInfo.m_initialState ) {

}

/*! \class Reference3d
 * Class for the GNDS **reference** distribution node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the Reference3d.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

Reference3d::Reference3d( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::reference3d, a_parent ),
        m_href( a_node.attribute_as_string( GIDI_hrefChars ) ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Reference3d::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;

    attributes += a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
    attributes += a_writeInfo.addAttribute( GIDI_hrefChars, href( ) );
    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}

/*! \class CoulombPlusNuclearElastic
 * Class for the GNDS **CoulombPlusNuclearElastic** distribution node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the CoulombPlusNuclearElastic.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

CoulombPlusNuclearElastic::CoulombPlusNuclearElastic( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::CoulombPlusNuclearElastic3d, a_parent ),
        m_href( a_node.attribute_as_string( GIDI_hrefChars ) ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void CoulombPlusNuclearElastic::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;

    attributes += a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
    attributes += a_writeInfo.addAttribute( GIDI_hrefChars, href( ) );
    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}

/*! \class LLNLLegendre
 * Class for the LLNL/GNDS **LLNLLegendre** distribution node.
 */

/* *********************************************************************************************************//**
 * This class is woefully inadequate but some form is needed by the method Product::isCompleteParticle.
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the LLNLLegendre.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

LLNLLegendre::LLNLLegendre( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::LLNLLegendre, a_parent ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void LLNLLegendre::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes;

    attributes += a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), attributes );
}

/*! \class Unspecified
 * Class for the GNDS **unspecified** node.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the Unspecified.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    The **m_distribution** member of GIDI::Product this distribution form belongs to.
 ***********************************************************************************************************/

Unspecified::Unspecified( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Distribution( a_node, a_setupInfo, FormType::unspecified, a_parent ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/
 
void Unspecified::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    toXMLNodeStarter( a_writeInfo, a_indent );
    a_writeInfo.addNodeEnder( moniker( ) );
}

}                     // End of namespace Distributions.

}                     // End of namespace GIDI.
