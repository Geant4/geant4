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

namespace Styles {

/*! \class Suite
 * This is essentially the GIDI::Suite class with the addition of the **findLabelInLineage** method.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Suite::Suite( ) :
        GIDI::Suite( GIDI_stylesChars ) {

}

/* *********************************************************************************************************//**
 * Searches the Suite *a_suite* for a form with label *a_label* or, if not found, recursively ascends the **derivedFrom** until
 * a derived form is found. The *this* instance must be an <**styles**> node so that the **derivedFrom**s can be ascended. 
 * If no form is found, an empty string is returned.
 *
 * @param a_suite       [in]    The Suite, typically a component, whose forms are searched for a form with label *a_label* or one of its **derivedFrom**.
 * @param a_label       [in]    The label of the form to start the search.
 * @return                      The label of the form found or an empty string if none is found.
 ***********************************************************************************************************/

std::string const *Suite::findLabelInLineage( GIDI::Suite const &a_suite, std::string const &a_label ) const {

    std::string const *label = &a_label;
    Suite::const_iterator iter = a_suite.find( a_label );

    while( true ) {
        if( iter != a_suite.end( ) ) return( label );

        Base const *form = get<Base>( *label );
        form = form->getDerivedStyle( );
        label = &form->keyValue( );
        if( *label == "" ) break;
        iter = a_suite.find( *label );        
    }

    return( label );
}

/*! \class Base
 * This is the virtual base class inherited by all **style** classes. It handles the *date* and **derivedFrom** members.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent      [in]    The parent GIDI::Suite.
 * @return
 ***********************************************************************************************************/

Base::Base( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent ) : 
        Form( a_node, a_setupInfo, FormType::style, a_parent ),
        m_date( a_node.attribute_as_string( GIDI_dateChars ) ),
        m_label( a_node.attribute_as_string( GIDI_labelChars ) ),
        m_derivedStyle( a_node.attribute_as_string( GIDI_derivedFromChars ) ) {

    if( a_node.child( GUPI_documentationChars ).empty( ) ) {
        m_documentation = nullptr;
    } else {
        m_documentation = new GUPI::Documentation( a_node.child( GIDI_documentationChars ) );
    }
}

Base::~Base( ) {

    delete m_documentation;

}

/* *********************************************************************************************************//**
 * Returns a pointer to the **derivedFrom** style of *this*.
 *
 * @return          Pointer to the **derivedFrom** style of *this*.
 ***********************************************************************************************************/

Base const *Base::getDerivedStyle( ) const {

    Form const *_form( sibling( m_derivedStyle ) );

    return( dynamic_cast<Base const *>( _form ) );
}

/* *********************************************************************************************************//**
 * Starting at *this*'s **derivedFrom** style, and ascending as needed, returns the **derivedFrom** style whose moniker is *a_moniker*.
 *
 * @param a_moniker         [in]    The moniker to search for.
 * @return                          The style whose moniker is *a_moniker*.
 ***********************************************************************************************************/

Base const *Base::getDerivedStyle( std::string const &a_moniker ) const {

    Form const *_form( sibling( m_derivedStyle ) );
    Base const *_style = dynamic_cast<Base const *>( _form );

    if( _style == nullptr ) return( _style );
    if( _style->moniker( ) != a_moniker ) _style = _style->getDerivedStyle( a_moniker );
    return( _style );
}

/* *********************************************************************************************************//**
 * Returns the base attributes for *this* as a *std::string* instance.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 *
 * @return                                      The base attributes as a XML attribute string.
 ***********************************************************************************************************/

std::string Base::baseXMLAttributes( GUPI::WriteInfo &a_writeInfo ) const {

    std::string attributes( a_writeInfo.addAttribute( GIDI_labelChars, label( ) ) );

    if( m_derivedStyle != "" ) attributes += a_writeInfo.addAttribute( GIDI_derivedFromChars, m_derivedStyle );
    attributes += a_writeInfo.addAttribute( GIDI_dateChars, m_date );

    return( attributes );
}

/*! \class Evaluated
 * This is the **evaluated** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent      [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Evaluated::Evaluated( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ),
        m_library( a_node.attribute_as_string( GIDI_libraryChars ) ),
        m_version( a_node.attribute_as_string( GIDI_versionChars ) ),
        m_temperature( a_node.child( GIDI_temperatureChars ), a_setupInfo ),
        m_projectileEnergyDomain( a_node.child( GIDI_projectileEnergyDomainChars ), a_setupInfo ) {

}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Evaluated::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes = baseXMLAttributes( a_writeInfo );

    attributes += a_writeInfo.addAttribute( GIDI_libraryChars, m_library );
    attributes += a_writeInfo.addAttribute( GIDI_versionChars, m_version );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    
    m_temperature.toXMLList( a_writeInfo, indent2 );
    m_projectileEnergyDomain.toXMLList( a_writeInfo, indent2 );
    
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class CrossSectionReconstructed
 * This is the **crossSectionReconstructed** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent      [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

CrossSectionReconstructed::CrossSectionReconstructed( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ),
        m_temperature( nullptr ) {

    HAPI::Node const temperatureNode = a_node.child( GIDI_temperatureChars );
    if( !temperatureNode.empty( ) ) {
        m_temperature = new PhysicalQuantity( temperatureNode, a_setupInfo );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

CrossSectionReconstructed::~CrossSectionReconstructed( ) {

    delete m_temperature;
}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/

PhysicalQuantity const &CrossSectionReconstructed::temperature( ) const {

    if( m_temperature != nullptr ) return( *m_temperature );

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void CrossSectionReconstructed::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), baseXMLAttributes( a_writeInfo ) );
    if( m_temperature != nullptr ) m_temperature->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class AngularDistributionReconstructed
 * This is the **angularDistributionReconstructed** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent      [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

AngularDistributionReconstructed::AngularDistributionReconstructed( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ),
        m_temperature( nullptr ) {

    HAPI::Node const temperatureNode = a_node.child( GIDI_temperatureChars );
    if( !temperatureNode.empty( ) ) {
        m_temperature = new PhysicalQuantity( temperatureNode, a_setupInfo );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

AngularDistributionReconstructed::~AngularDistributionReconstructed( ) {

    delete m_temperature;
}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/

PhysicalQuantity const &AngularDistributionReconstructed::temperature( ) const {

    if( m_temperature != nullptr ) return( *m_temperature );

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void AngularDistributionReconstructed::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), baseXMLAttributes( a_writeInfo ) );
    if( m_temperature != nullptr ) m_temperature->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class CoulombPlusNuclearElasticMuCutoff
 * This is the **CoulombPlusNuclearElasticMuCutoff** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent      [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

CoulombPlusNuclearElasticMuCutoff::CoulombPlusNuclearElasticMuCutoff( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ),
        m_muCutoff( a_node.attribute_as_double( GIDI_muCutoffChars ) ) {

}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/

PhysicalQuantity const &CoulombPlusNuclearElasticMuCutoff::temperature( ) const {

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void CoulombPlusNuclearElasticMuCutoff::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes = baseXMLAttributes( a_writeInfo );

    attributes += a_writeInfo.addAttribute( GIDI_muCutoffChars, LUPI::Misc::doubleToShortestString( m_muCutoff ) );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class Realization
 * This is the GNDS **Realization** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent      [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Realization::Realization( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ) {

}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/

PhysicalQuantity const &Realization::temperature( ) const {

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Realization::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string attributes = baseXMLAttributes( a_writeInfo );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class AverageProductData
 * This is the **averageProductData** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent      [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

AverageProductData::AverageProductData( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ),
        m_temperature( nullptr ) {

    HAPI::Node const temperatureNode = a_node.child( GIDI_temperatureChars );
    if( !temperatureNode.empty( ) ) {
        m_temperature = new PhysicalQuantity( temperatureNode, a_setupInfo );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

AverageProductData::~AverageProductData( ) {

    delete m_temperature;
}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/

PhysicalQuantity const &AverageProductData::temperature( ) const {

    if( m_temperature != nullptr ) return( *m_temperature );

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void AverageProductData::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), baseXMLAttributes( a_writeInfo ) );
    if( m_temperature != nullptr ) m_temperature->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class Heated
 * This is the **heated** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent      [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Heated::Heated( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ),
        m_temperature( a_node.child( GIDI_temperatureChars ), a_setupInfo ) {
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Heated::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), baseXMLAttributes( a_writeInfo ) );
    m_temperature.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class MonteCarlo_cdf
 * This is the **MonteCarlo_cdf** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent      [in]    The parent GIDI::Suite.
 * @return
 ***********************************************************************************************************/

MonteCarlo_cdf::MonteCarlo_cdf( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ) {

}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/

PhysicalQuantity const &MonteCarlo_cdf::temperature( ) const {

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void MonteCarlo_cdf::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    a_writeInfo.addNodeStarter( a_indent, moniker( ), baseXMLAttributes( a_writeInfo ) );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class MultiGroup
 * This is the **multiGroup** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

MultiGroup::MultiGroup( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops,
		        PoPI::Database const &a_internalPoPs, GIDI::Suite *a_parent ) : 
        Base( a_node, a_setupInfo, a_parent ),
        m_maximumLegendreOrder( a_node.attribute_as_int( GIDI_lMaxChars ) ),
        m_transportables( a_construction, GIDI_transportablesChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseTransportablesSuite, nullptr ) {

    m_transportables.setAncestor( this );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

MultiGroup::~MultiGroup( ) {

}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/
 
PhysicalQuantity const &MultiGroup::temperature( ) const {

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Returns the multi-group boundaries for the product with index *a_productID*.
 *
 * @param a_productID           [in]    Particle id for the requested product.
 * @return                              The multi-group boundaries.
 ***********************************************************************************************************/

std::vector<double> MultiGroup::groupBoundaries( std::string const &a_productID ) const {

    for( std::size_t index = 0; index < m_transportables.size( ); ++index ) {
        Transportable const &transportable1 = *m_transportables.get<Transportable>( index );

        if( transportable1.pid( ) == a_productID ) {
            return( transportable1.groupBoundaries( ) );
        }
    }
    throw Exception( "MultiGroup::groupBoundaries: product index not found" );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void MultiGroup::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes = baseXMLAttributes( a_writeInfo );

    attributes += a_writeInfo.addAttribute( GIDI_lMaxChars, intToString( m_maximumLegendreOrder ) );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    m_transportables.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class HeatedMultiGroup
 * This is the **neatedMultiGroup** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

HeatedMultiGroup::HeatedMultiGroup( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, 
                PoPI::Database const &a_pops, GIDI::Suite *a_parent ) : 
        Base( a_node, a_setupInfo, a_parent ),
        m_transportables( a_construction, GIDI_transportablesChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_pops, parseTransportablesSuite, nullptr ),
        m_flux( a_construction, a_node.child( GIDI_fluxNodeChars ), a_setupInfo ),
        m_inverseSpeed( a_construction, a_node.child( GIDI_inverseSpeedChars ).child( GIDI_gridded1dChars ), a_setupInfo, nullptr ),
        m_parameters( a_node.attribute_as_string( GIDI_parametersChars ) ) {

    m_transportables.setAncestor( this );
    m_flux.setAncestor( this );
    m_inverseSpeed.setAncestor( this );

    if( m_transportables.size( ) == 0 ) {
        GIDI::Suite const *transportables1 = nullptr;
        if( a_setupInfo.m_multiGroup != nullptr ) {
            transportables1 = &a_setupInfo.m_multiGroup->transportables( ); }
        else if( a_setupInfo.m_heatedMultiGroup != nullptr ) {
            transportables1 = &a_setupInfo.m_heatedMultiGroup->transportables( );
        }
        if( transportables1 != nullptr ) {
            for( std::size_t index = 0; index < transportables1->size( ); ++index ) {
                Transportable const &transportable = *transportables1->get<Transportable>( index );

                m_transportables.add( new Transportable( transportable ) );
            }
            
        } }
    else if( a_setupInfo.m_heatedMultiGroup == nullptr ) {
        a_setupInfo.m_heatedMultiGroup = this;
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

HeatedMultiGroup::~HeatedMultiGroup( ) {

}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/

PhysicalQuantity const &HeatedMultiGroup::temperature( ) const {

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Returns the **Transportable** instance for the particle with id *a_ID* used for processing this **HeatedMultiGroup**.
 *
 * @param a_ID                  [in]    Particle id for the requested product.
 * @return                              The multi-group boundaries.
 ***********************************************************************************************************/

Transportable const &HeatedMultiGroup::transportable( std::string const &a_ID ) const {

    return( *m_transportables.get<Transportable>( a_ID ) );
}


/* *********************************************************************************************************//**
 * Returns the multi-group boundaries for the particle with id *a_ID* used for processing this **HeatedMultiGroup**.
 *
 * @param a_ID                  [in]    Particle id for the requested product.
 * @return                              The multi-group boundaries.
 ***********************************************************************************************************/

std::vector<double> HeatedMultiGroup::groupBoundaries( std::string const &a_ID ) const {

    Transportable const &transportable1 = transportable( a_ID );

    return( transportable1.groupBoundaries( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void HeatedMultiGroup::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string indent3 = a_writeInfo.incrementalIndent( indent2 );
    std::string attributes = baseXMLAttributes( a_writeInfo );
    
    attributes += a_writeInfo.addAttribute( GIDI_parametersChars, m_parameters );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    m_flux.toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeStarter( indent2, GIDI_inverseSpeedChars, "" );
    m_inverseSpeed.toXMLList_func( a_writeInfo, indent3, false, false );
    a_writeInfo.addNodeEnder( GIDI_inverseSpeedChars );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class SnElasticUpScatter
 * This is the **SnElasticUpScatter** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node        [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops        [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_parent      [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

SnElasticUpScatter::SnElasticUpScatter( HAPI::Node const &a_node, SetupInfo &a_setupInfo, LUPI_maybeUnused PoPI::Database const &a_pops, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ),
        m_upperCalculatedGroup( a_node.attribute_as_int( GIDI_upperCalculatedGroupChars ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

SnElasticUpScatter::~SnElasticUpScatter( ) {

}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/

PhysicalQuantity const &SnElasticUpScatter::temperature( ) const {

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void SnElasticUpScatter::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string attributes = baseXMLAttributes( a_writeInfo );
    
    attributes += a_writeInfo.addAttribute( GIDI_upperCalculatedGroupChars, intToString( m_upperCalculatedGroup ) );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class GriddedCrossSection
 * This is the **griddedCrossSection** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

GriddedCrossSection::GriddedCrossSection( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, LUPI_maybeUnused PoPI::Database const &a_pops, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ),
        m_grid( a_node.child( GIDI_gridChars ), a_setupInfo, a_construction.useSystem_strtod( ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

GriddedCrossSection::~GriddedCrossSection( ) {

}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/

PhysicalQuantity const &GriddedCrossSection::temperature( ) const {

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void GriddedCrossSection::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    
    a_writeInfo.addNodeStarter( a_indent, moniker( ), baseXMLAttributes( a_writeInfo ) );
    m_grid.toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/*! \class URR_probabilityTables
 * This is the **URR_probabilityTables** style class.
 */

/* *********************************************************************************************************//**
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    A PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

URR_probabilityTables::URR_probabilityTables( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, LUPI_maybeUnused PoPI::Database const &a_pops, GIDI::Suite *a_parent ) :
        Base( a_node, a_setupInfo, a_parent ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

URR_probabilityTables::~URR_probabilityTables( ) {

}

/* *********************************************************************************************************//**
 * Ascends the **derivedFrom** styles until a temperature is found.
 *
 * @return          Returns the temperature associated with this style.
 ***********************************************************************************************************/

PhysicalQuantity const &URR_probabilityTables::temperature( ) const {

    Base const *style = getDerivedStyle( );

    if( style == nullptr ) throw Exception( "No style with temperature." );
    return( style->temperature( ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void URR_probabilityTables::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    a_writeInfo.addNodeStarterEnder( a_indent, moniker( ), baseXMLAttributes( a_writeInfo ) );
}

/*! \class TemperatureInfo
 * This class stores the labels for a given temperature for the **heatedCrossSection**, **griddedCrossSection**, **heatedMultiGroup** and
 * **SnElasticUpScatter** styles. If no style of a given process (e.g., **heatedCrossSection**) type exists, its label is an empty string.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

TemperatureInfo::TemperatureInfo( ) :
        m_temperature( -1.0, "K" ),
        m_heatedCrossSection( "" ),
        m_griddedCrossSection( "" ),
        m_URR_probabilityTables( "" ),
        m_heatedMultiGroup( "" ),
        m_SnElasticUpScatter( "" ) {

}

/* *********************************************************************************************************//**
 *
 * @param a_temperature             [in]    The temperature.
 * @param a_heatedCrossSection      [in]    The label for the **heatedCrossSection** style.
 * @param a_griddedCrossSection     [in]    The label for the **griddedCrossSection** style.
 * @param a_heatedMultiGroup        [in]    The label for the **heatedMultiGroup** style.
 * @param a_URR_probabilityTables   [in]    The label for the **URR_probabilityTables** style.
 * @param a_SnElasticUpScatter      [in]    The label for the **SnElasticUpScatter** style.
 ***********************************************************************************************************/

TemperatureInfo::TemperatureInfo( PhysicalQuantity const &a_temperature, std::string const &a_heatedCrossSection, std::string const &a_griddedCrossSection,
                std::string const &a_URR_probabilityTables, std::string const &a_heatedMultiGroup, std::string const &a_SnElasticUpScatter ) :
        m_temperature( a_temperature ),
        m_heatedCrossSection( a_heatedCrossSection ),
        m_griddedCrossSection( a_griddedCrossSection ),
        m_URR_probabilityTables( a_URR_probabilityTables ),
        m_heatedMultiGroup( a_heatedMultiGroup ),
        m_SnElasticUpScatter( a_SnElasticUpScatter ) {

}

/* *********************************************************************************************************//**
 * Prints information about *this* to std::cout.
 ***********************************************************************************************************/

void TemperatureInfo::print( ) const {

    std::cout << "temperature = " << m_temperature.value( ) << " " << m_temperature.unit( ) << " heatedCrossSection = '" << m_heatedCrossSection
            << "' griddedCrossSection = '" << m_griddedCrossSection << "' URR_probabilityTables = '" << m_URR_probabilityTables 
            << "' heatedMultiGroup = '" << m_heatedMultiGroup << "' SnElasticUpScatter = '" << m_SnElasticUpScatter << std::endl;
}

}

}
