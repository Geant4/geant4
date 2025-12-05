/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <cmath>

#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

/* *********************************************************************************************************//**
 * Parses a <**reaction**> node.
 ***********************************************************************************************************/

Reaction::Reaction( int a_ENDF_MT, std::string const &a_fissionGenre ) :
        Form( FormType::reaction ),
        m_reactionIndex( 0 ),
        m_active( true ),
        m_ENDF_MT( a_ENDF_MT ),
        m_fissionGenre( a_fissionGenre ),
        m_QThreshold( 0.0 ),
        m_crossSectionThreshold( 0.0 ),
        m_twoBodyThreshold( 0.0 ),
        m_isPairProduction( false ),
        m_isPhotoAtomicIncoherentScattering( false ),
        m_RutherfordScatteringPresent( false ),
        m_onlyRutherfordScatteringPresent( false ),
        m_nuclearPlusInterferencePresent( false ),
        m_decayPositronium( false ),
        m_doubleDifferentialCrossSection( GIDI_doubleDifferentialCrossSectionChars, GIDI_labelChars ),
        m_crossSection( GIDI_crossSectionChars, GIDI_labelChars ),
        m_availableEnergy( GIDI_availableEnergyChars, GIDI_labelChars ),
        m_availableMomentum( GIDI_availableMomentumChars, GIDI_labelChars ),
        m_outputChannel( ) {

    setMoniker( GIDI_reactionChars );

    ENDL_CFromENDF_MT( m_ENDF_MT, &m_ENDL_C, &m_ENDL_S );
}

/* *********************************************************************************************************//**
 * Parses a <**reaction**> node.
 *
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The reaction HAPI::Node to be parsed and used to construct the reaction.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops            [in]    The *external* PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs    [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                  This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_protare         [in]    The GIDI::Protare this reaction belongs to.
 * @param a_styles          [in]    The <**styles**> node under the <**reactionSuite**> node.
 ***********************************************************************************************************/

Reaction::Reaction( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, 
                PoPI::Database const &a_internalPoPs, Protare const &a_protare, Styles::Suite const *a_styles ) :
        Form( a_node, a_setupInfo, FormType::reaction ),
        m_active( true ),
        m_ENDF_MT( a_node.attribute_as_int( GIDI_ENDF_MT_Chars ) ),
        m_ENDL_C( 0 ),
        m_ENDL_S( 0 ),
        m_fissionGenre( a_node.attribute_as_string( GIDI_fissionGenreChars ) ),
        m_QThreshold( 0.0 ),
        m_crossSectionThreshold( 0.0 ),
        m_twoBodyThreshold( 0.0 ),
        m_isPairProduction( false ),
        m_isPhotoAtomicIncoherentScattering( false ),
        m_RutherfordScatteringPresent( false ),
        m_onlyRutherfordScatteringPresent( false ),
        m_nuclearPlusInterferencePresent( false ),
        m_decayPositronium( false ),
        m_doubleDifferentialCrossSection( a_construction, GIDI_doubleDifferentialCrossSectionChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, 
                a_internalPoPs, parseDoubleDifferentialCrossSectionSuite, a_styles ),
        m_crossSection( a_construction, GIDI_crossSectionChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseCrossSectionSuite, a_styles ),
        m_availableEnergy( a_construction, GIDI_availableEnergyChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseAvailableSuite, a_styles ),
        m_availableMomentum( a_construction, GIDI_availableMomentumChars, GIDI_labelChars, a_node, a_setupInfo, a_pops, a_internalPoPs, parseAvailableSuite, a_styles ),
        m_outputChannel( nullptr ) {

    m_isPairProduction = label( ).find( "pair production" ) != std::string::npos;
    m_decayPositronium = m_isPairProduction && a_construction.decayPositronium( );
    m_isPhotoAtomicIncoherentScattering = false;
    if( m_doubleDifferentialCrossSection.size( ) > 0 )
        m_isPhotoAtomicIncoherentScattering = m_doubleDifferentialCrossSection.get<Form>( 0 )->type( ) == FormType::incoherentPhotonScattering;

    m_doubleDifferentialCrossSection.setAncestor( this );
    m_crossSection.setAncestor( this );
    m_availableEnergy.setAncestor( this );
    m_availableMomentum.setAncestor( this );

    a_setupInfo.m_outputChannelLevel = 0;
    m_outputChannel = new OutputChannel( a_construction, a_node.child( GIDI_outputChannelChars ), a_setupInfo, a_pops, 
            a_internalPoPs, a_styles, hasFission( ), true );
    m_outputChannel->setAncestor( this );

    HAPI::Node const CoulombPlusNuclearElastic = a_node.child( GIDI_doubleDifferentialCrossSectionChars ).first_child( );
    if( CoulombPlusNuclearElastic.name( ) == GIDI_CoulombPlusNuclearElasticChars ) {    // Check RutherfordScattering.
        m_RutherfordScatteringPresent = true;
        m_onlyRutherfordScatteringPresent = true;
        for( HAPI::Node child = CoulombPlusNuclearElastic.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
            if( child.name( ) != GIDI_RutherfordScatteringChars ) m_onlyRutherfordScatteringPresent = false;
            if( child.name( ) == GIDI_nuclearPlusInterferenceChars ) m_nuclearPlusInterferencePresent = true;
        }
    }

    ENDL_CFromENDF_MT( m_ENDF_MT, &m_ENDL_C, &m_ENDL_S );
    if( a_setupInfo.m_isENDL_C_9 ) m_ENDL_C = 9;
    if( m_ENDL_C == 20 ) {
        Product *product = m_outputChannel->products( ).get<Product>( 0 );
        if( product->particle( ).ID( ) == "H1" ) m_ENDL_C = 21;
    }
    else if( m_ENDL_C == 22 ) {
        Product *product = m_outputChannel->products( ).get<Product>( 0 );
        if( product->particle( ).ID( ) == "H2" ) m_ENDL_C = 35;
    }

    if( ( a_construction.parseMode( ) != Construction::ParseMode::outline ) && ( a_construction.parseMode( ) != Construction::ParseMode::readOnly ) ) {
        double _Q = 0.0;
        Form const &QForm = *(m_outputChannel->Q( ).get<Form>( 0 ));
        switch( QForm.type( ) ) {
        case FormType::constant1d : {
            Functions::Constant1d const &constant1d = static_cast<Functions::Constant1d const &>( QForm );
            _Q = constant1d.value( );
            break; }
        case FormType::XYs1d : {
            Functions::XYs1d const &xys1d = static_cast<Functions::XYs1d const &>( QForm );
            _Q = xys1d.evaluate( 0.0 );
            break; }
        case FormType::gridded1d : {        // Should be a special reaction (e.g., summed multi-group data) and can be ignored.
            _Q = 0.0;
            break; }
        default :
            throw Exception( "Reaction::Reaction: unsupported Q form " + QForm.label( ) );
        }
        m_twoBodyThreshold = -_Q * a_protare.thresholdFactor( );
        _Q *= -1;
        if( _Q <= 0.0 ) _Q = 0.0;
        m_QThreshold = a_protare.thresholdFactor( ) * _Q;

        if( _Q > 0.0 ) {               // Try to do a better job determining m_crossSectionThreshold.
            std::vector<Suite::const_iterator> monikers = a_protare.styles( ).findAllOfMoniker( GIDI_griddedCrossSectionStyleChars );
            if( ( a_construction.parseMode( ) != Construction::ParseMode::multiGroupOnly ) && ( monikers.size( ) > 0 ) ) {
                Styles::GriddedCrossSection const &griddedCrossSection = static_cast<Styles::GriddedCrossSection const &>( **monikers[0] );
                Grid grid = griddedCrossSection.grid( );

                if( m_crossSection.has( griddedCrossSection.label( ) ) ) {
                    Functions::Ys1d const &ys1d = static_cast<Functions::Ys1d const &>( *m_crossSection.get<Functions::Ys1d>( griddedCrossSection.label( ) ) );
                    m_crossSectionThreshold = grid[ys1d.start( )];
                }
            }
            if( m_crossSectionThreshold == 0.0 ) {      // Should also check 'evaluate' style before using m_QThreshold as a default.
                m_crossSectionThreshold = m_QThreshold;
            }
        }
        if( m_twoBodyThreshold > 0.0 ) m_twoBodyThreshold = m_crossSectionThreshold;
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Reaction::~Reaction( ) {

    if( m_outputChannel != nullptr ) delete m_outputChannel;
}

/* *********************************************************************************************************//**
 * Fills in a std::set with a unique list of all product indices produced by this reaction.
 * If a_transportablesOnly is true, only transportable product indices are return.
 *
 * @param a_ids                     [out]   The unique list of product indices.
 * @param a_particles               [in]    The list of particles to be transported.
 * @param a_transportablesOnly      [in]    If true, only transportable product indices are added in the list.
 ***********************************************************************************************************/

void Reaction::productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const {

    if( m_decayPositronium ) {
        std::string electronAnti( PoPI::IDs::electron + PoPI::IDs::anti );
        std::set<std::string> ids;
        m_outputChannel->productIDs( ids, a_particles, a_transportablesOnly );
        for( auto iterId = ids.begin( ); iterId != ids.end( ); ++iterId ) {
            if( ( *iterId == PoPI::IDs::electron ) || ( *iterId == electronAnti ) ) continue;
            a_ids.insert( *iterId );
        }
        a_ids.insert( PoPI::IDs::photon ); }
    else {
        m_outputChannel->productIDs( a_ids, a_particles, a_transportablesOnly );
    }
}

/* *********************************************************************************************************//**
 * Determines the maximum Legredre order present in the multi-group transfer matrix for a give product for a give label. Inspects all
 * products produced by this reaction.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    Particle id of the requested product.
 *
 * @return                              The maximum Legredre order. If no transfer matrix data are present for the requested product, -1 is returned.
 ***********************************************************************************************************/

int Reaction::maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    if( m_isPairProduction && ( a_productID == PoPI::IDs::photon ) ) return( 0 );
    return( m_outputChannel->maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID ) );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total multiplicity for the requested label for the requested product. This is a cross section weighted multiplicity.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    Particle id for the requested product.
 *
 * @return                              The requested multi-group multiplicity as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Reaction::multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                    Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    if( m_decayPositronium ) {
        if( a_productID == PoPI::IDs::photon ) vector += multiGroupCrossSection( a_smr, a_settings, a_temperatureInfo ) * 2; }
    else {
        vector += m_outputChannel->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns true if at least one output channel contains a fission channel.
 *
 * @return  true if at least one output channel is a fission channel.
 ***********************************************************************************************************/

bool Reaction::hasFission( ) const {

    if( m_ENDF_MT == 18 ) return( true );
    if( m_ENDF_MT == 19 ) return( true );
    if( m_ENDF_MT == 20 ) return( true );
    if( m_ENDF_MT == 21 ) return( true );
    if( m_ENDF_MT == 38 ) return( true );

    return( false );
}

/* *********************************************************************************************************//**
 * Sets *this* reaction's output channel to **a_outputChannel**.
 *
 * @param a_outputChannel   [in]    The output channel to make *this* reaction output channel.
 ***********************************************************************************************************/

void Reaction::setOutputChannel( OutputChannel *a_outputChannel ) {

    m_outputChannel = a_outputChannel;
    m_outputChannel->setAncestor( this );
}

/* *********************************************************************************************************//**
 * Only for internal use. Called by ProtareTNSL instance to zero the lower energy multi-group data covered by the TNSL ProtareSingle.
 *
 * @param a_maximumTNSL_MultiGroupIndex     [in]    A map that contains labels for heated multi-group data and the last valid group boundary
 *                                                  for the TNSL data for that boundary.
 ***********************************************************************************************************/

void Reaction::modifiedMultiGroupElasticForTNSL( std::map<std::string,std::size_t> const &a_maximumTNSL_MultiGroupIndex ) {

    m_crossSection.modifiedMultiGroupElasticForTNSL( a_maximumTNSL_MultiGroupIndex );
    m_availableEnergy.modifiedMultiGroupElasticForTNSL( a_maximumTNSL_MultiGroupIndex );
    m_availableMomentum.modifiedMultiGroupElasticForTNSL( a_maximumTNSL_MultiGroupIndex );
    m_outputChannel->modifiedMultiGroupElasticForTNSL( a_maximumTNSL_MultiGroupIndex );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 *
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry *Reaction::findInAncestry3( std::string const &a_item ) {

    if( a_item == GIDI_doubleDifferentialCrossSectionChars ) return( &m_doubleDifferentialCrossSection );
    if( a_item == GIDI_crossSectionChars ) return( &m_crossSection );
    if( a_item == GIDI_availableEnergyChars ) return( &m_availableEnergy );
    if( a_item == GIDI_availableMomentumChars ) return( &m_availableMomentum );
    if( a_item == GIDI_outputChannelChars ) return( m_outputChannel );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 *
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry const *Reaction::findInAncestry3( std::string const &a_item ) const {

    if( a_item == GIDI_doubleDifferentialCrossSectionChars ) return( &m_doubleDifferentialCrossSection );
    if( a_item == GIDI_crossSectionChars ) return( &m_crossSection );
    if( a_item == GIDI_availableEnergyChars ) return( &m_availableEnergy );
    if( a_item == GIDI_availableMomentumChars ) return( &m_availableMomentum );
    if( a_item == GIDI_outputChannelChars ) return( m_outputChannel );

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns **true** if all outgoing particles (i.e., products) are specifed in *a_particles*. That is, the user
 * will be tracking all products of *this* reaction.
 *
 * @param a_particles           [in]    The list of particles to be transported.
 *
 * @return                              bool.
 ***********************************************************************************************************/

bool Reaction::areAllProductsTracked( Transporting::Particles const &a_particles ) const {

    if( hasFission( ) || m_isPhotoAtomicIncoherentScattering ) return( false );

    return( m_outputChannel->areAllProductsTracked( a_particles ) );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, cross section for the requested label the reaction.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_label               [in]    If not an empty string, this is used as the label for the form to return and the *a_temperatureInfo* labels are ignored.
 *
 * @return                              The requested multi-group cross section as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Reaction::multiGroupCrossSection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_label ) const {

    Vector vector( 0 );

    Functions::Gridded1d const *form = dynamic_cast<Functions::Gridded1d const*>( 
            a_settings.form( a_smr, m_crossSection, a_temperatureInfo, "cross section", a_label ) );
    if( form != nullptr ) vector = form->data( );

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total Q for the requested label. This is a cross section weighted Q.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_final               [in]    If false, only the Q for the primary reactions are return, otherwise, the Q for the final reactions.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the Q.
 *
 * @return                          The requested multi-group Q as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Reaction::multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings,
                Styles::TemperatureInfo const &a_temperatureInfo, bool a_final ) const {

    if( m_decayPositronium ) return( Vector { 0 } );          // Special case, returns Q with all 0.0s.

    return( m_outputChannel->multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final ) );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, product matrix for the requested label for the requested product index for the requested Legendre order.
 * If no data are found, an empty GIDI::Matrix is returned.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label and if delayed neutrons should be included.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_productID           [in]    Particle id for the requested product.
 * @param a_order               [in]    Requested product matrix, Legendre order.
 *
 * @return                              The requested multi-group product matrix as a GIDI::Matrix.
 ***********************************************************************************************************/

Matrix Reaction::multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::string const &a_productID, 
                std::size_t a_order ) const {

    Matrix matrix( 0, 0 );

    if( m_decayPositronium ) {
        if( a_productID == PoPI::IDs::photon ) {
            if( a_order == 0 ) {
                Vector productionCrossSection = multiGroupCrossSection( a_smr, a_settings, a_temperatureInfo ) * 2;
                std::map<std::string, GIDI::Transporting::Particle> const &particles = a_particles.particles( );
                std::map<std::string, GIDI::Transporting::Particle>::const_iterator particle = particles.find( PoPI::IDs::photon );
                GIDI::Transporting::MultiGroup const &multiGroup = particle->second.fineMultiGroup( );
                int multiGroupIndexFromEnergy = multiGroup.multiGroupIndexFromEnergy( PoPI_electronMass_MeV_c2, true );
                Matrix matrix2( productionCrossSection.size( ), productionCrossSection.size( ) );

                for( std::size_t i1 = 0; i1 < productionCrossSection.size( ); ++i1 ) {
                    matrix2.set( i1, static_cast<std::size_t>( multiGroupIndexFromEnergy ), productionCrossSection[i1] );
                }
                matrix += matrix2;
            }
        } }
    else {
        matrix += m_outputChannel->multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order );
    }

    return( matrix );
}

/* *********************************************************************************************************//**
 * Like Reaction::multiGroupProductMatrix, but only returns the fission neutron, transfer matrix.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label and if delayed neutrons should be included.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_order               [in]    Requested product matrix, Legendre order.
 *
 * @return                              The requested multi-group neutron fission matrix as a GIDI::Matrix.
 ***********************************************************************************************************/

Matrix Reaction::multiGroupFissionMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::size_t a_order ) const {

    Matrix matrix( 0, 0 );

    if( hasFission( ) ) matrix += multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, PoPI::IDs::neutron, a_order );
    return( matrix );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total available energy for the requested label. This is a cross section weighted available energy.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 *
 * @return                              The requested multi-group available energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Reaction::multiGroupAvailableEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo ) const {

    Vector vector( 0 );

    Functions::Gridded1d const *form = dynamic_cast<Functions::Gridded1d const*>( a_settings.form( a_smr, m_availableEnergy, a_temperatureInfo, "available energy" ) );
    if( form != nullptr ) vector = form->data( );

    if( m_decayPositronium ) 
        vector -= m_outputChannel->multiGroupQ( a_smr, a_settings, a_temperatureInfo, false );

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total average energy for the requested label for the requested product. This is a cross section weighted average energy
 * summed over all products for this reaction.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    Particle id for the requested product.
 *
 * @return                              The requested multi-group average energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Reaction::multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    if( m_decayPositronium ) {
        if( a_productID == PoPI::IDs::photon ) vector += multiGroupCrossSection( a_smr, a_settings, a_temperatureInfo ) * 2.0 * PoPI_electronMass_MeV_c2; }
    else {
        vector += m_outputChannel->multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total deposition energy for the requested label for *this* reaction. This is a cross section weighted 
 * deposition energy. The deposition energy is calculated by subtracting the average energy from each transportable particle
 * from the available energy. The list of transportable particles is specified via the list of particle specified in the *a_settings* argument.
 * This method does not include any photon deposition energy for *this* reaction that is in the **GNDS** orphanProducts node.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label and the products that are transported.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 *
 * @return                              The requested multi-group deposition energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Reaction::multiGroupDepositionEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles ) const {

    std::map<std::string, Transporting::Particle> const &products = a_particles.particles( );
    Vector vector;

    if( moniker( ) != GIDI_orphanProductChars ) {
        if( ( a_settings.zeroDepositionIfAllProductsTracked( ) ) && areAllProductsTracked( a_particles ) && !m_isPairProduction ) return( vector );

        vector = multiGroupAvailableEnergy( a_smr, a_settings, a_temperatureInfo );
    }

    Vector availableEnergy( vector );

    for( std::map<std::string, Transporting::Particle>::const_iterator iter = products.begin( ); iter != products.end( ); ++iter ) {
        vector -= multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, iter->first );
    }

    for( std::size_t index = 0; index < availableEnergy.size( ); ++index ) {        // Check for values that should probably be 0.0.
        if( std::fabs( vector[index] ) < std::fabs( 1e-14 * availableEnergy[index] ) ) vector[index] = 0.0;
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total available momentum for the requested label. This is a cross section weighted available momentum.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 *
 * @return                              The requested multi-group available momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Reaction::multiGroupAvailableMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo ) const {

    Vector vector( 0 );

    Functions::Gridded1d const *form = dynamic_cast<Functions::Gridded1d const*>( a_settings.form( a_smr, m_availableMomentum, a_temperatureInfo, "available momentum" ) );
    if( form != nullptr ) vector = form->data( );

    return( vector );
}
/* *********************************************************************************************************//**
 * Returns the multi-group, total average momentum for the requested label for the requested product. This is a cross section weighted average momentum.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    Particle id for the requested product.
 *
 * @return                              The requested multi-group average momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Reaction::multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    Vector vector( 0 );

    if( !m_isPairProduction ) vector += m_outputChannel->multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID );

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total deposition momentum for the requested label for *this* reaction. This is a cross section 
 * weighted deposition momentum. The deposition momentum is calculated by subtracting the average momentum from each transportable particle
 * from the available momentum. The list of transportable particles is specified via the list of particle specified in the *a_settings* argument.
 * This method does not include any photon deposition momentum for *this* reaction that is in the **GNDS** orphanProducts node.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 *
 * @return                              The requested multi-group deposition momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector Reaction::multiGroupDepositionMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles ) const {

    std::map<std::string, Transporting::Particle> const &products = a_particles.particles( );
    Vector vector;

    if( moniker( ) != GIDI_orphanProductChars ) {
        vector = multiGroupAvailableMomentum( a_smr, a_settings, a_temperatureInfo );
    }

    for( std::map<std::string, Transporting::Particle>::const_iterator iter = products.begin( ); iter != products.end( ); ++iter ) {
        vector -= multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, iter->first );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, gain for the requested particle and label. This is a cross section weighted gain summed over all reactions.
 * If *a_productID* and *a_projectileID* are the same, then the multi-group cross section is subtracted for the returned value to indicate
 * that the *a_projectileID* as been absorted.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    The particle PoPs' id for the whose gain is to be calculated.
 * @param a_projectileID        [in]    The particle PoPs' id for the projectile.
 *
 * @return                              The requested multi-group gain as a **GIDI::Vector**.
 ***********************************************************************************************************/

Vector Reaction::multiGroupGain( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID, std::string const &a_projectileID  ) const {

    Vector vector( multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID ) );

    if( PoPI::compareSpecialParticleIDs( a_productID, a_projectileID ) ) vector -= multiGroupCrossSection( a_smr, a_settings, a_temperatureInfo );

    return( vector );
}

/* *********************************************************************************************************//**
 * Appends a DelayedNeutronProduct instance for each delayed neutron in *m_delayedNeutrons*.
 *
 * @param       a_delayedNeutronProducts    [in/out]    The list to append the delayed neutrons to.
 ***********************************************************************************************************/

void Reaction::delayedNeutronProducts( DelayedNeutronProducts &a_delayedNeutronProducts ) const {

    if( m_outputChannel != nullptr ) m_outputChannel->delayedNeutronProducts( a_delayedNeutronProducts );
}

/* *********************************************************************************************************//**
 * If *this* has an output channel, this output channel's **incompleteParticles* is called.
 *
 * @param a_settings                    [in]    Specifies the requested label.
 * @param a_incompleteParticles         [out]   The list of particles whose **completeParticle** method returns *false*.
 ***********************************************************************************************************/

void Reaction::incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const {

    if( m_outputChannel != nullptr ) m_outputChannel->incompleteParticles( a_settings, a_incompleteParticles );
    if( m_isPairProduction ) a_incompleteParticles.erase( PoPI::IDs::photon );          // Kludge for old processed files.
}

/* *********************************************************************************************************//**
 * Returns, via arguments, the average energy and momentum, and gain for product with particle id *a_particleID*.
 *
 * @param a_settings                    [in]    Specifies the requested label.
 * @param a_particleID                  [in]    The particle id of the product.
 * @param a_energy                      [in]    The energy of the projectile.
 * @param a_productEnergy               [in]    The average energy of the product.
 * @param a_productMomentum             [in]    The average momentum of the product.
 * @param a_productGain                 [in]    The gain of the product.
 * @param a_ignoreIncompleteParticles   [in]    If *true*, incomplete particles are ignore, otherwise a *throw* is executed.
 ***********************************************************************************************************/

void Reaction::continuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, double a_energy, double &a_productEnergy, double &a_productMomentum, 
                double &a_productGain, bool a_ignoreIncompleteParticles ) const {

    a_productEnergy = 0.0;
    a_productMomentum = 0.0;
    a_productGain = 0.0;

//    if( ENDF_MT( ) == 516 ) return;             // FIXME, may be something wrong with the way FUDGE converts ENDF to GNDS.

    if( m_outputChannel != nullptr )
        m_outputChannel->continuousEnergyProductData( a_settings, a_particleID, a_energy, a_productEnergy, a_productMomentum, a_productGain, 
                a_ignoreIncompleteParticles );
}

/* *********************************************************************************************************//**
 * Modifies the average product energies, momenta and gains for product with particle id *a_particleID*.
 *
 * @param a_settings                    [in]    Specifies user options.
 * @param a_particleID                  [in]    The particle id of the product.
 * @param a_energies                    [in]    The vector of energies to map the data to.
 * @param a_offset                      [in]    The index of the first energy whose data are to be added to the vectors.
 * @param a_productEnergies             [out]   The vector of average energies of the product.
 * @param a_productMomenta              [out]   The vector of average momenta of the product.
 * @param a_productGains                [out]   The vector of gain of the product.
 * @param a_ignoreIncompleteParticles   [in]    If *true*, incomplete particles are ignore, otherwise a *throw* is executed.
 ***********************************************************************************************************/

void Reaction::mapContinuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, 
                std::vector<double> const &a_energies, std::size_t a_offset, std::vector<double> &a_productEnergies, std::vector<double> &a_productMomenta, 
                std::vector<double> &a_productGains, bool a_ignoreIncompleteParticles ) const {

//    if( ENDF_MT( ) == 516 ) return;             // FIXME, may be something wrong with the way FUDGE converts ENDF to GNDS.

    if( m_outputChannel != nullptr )
        m_outputChannel->mapContinuousEnergyProductData( a_settings, a_particleID, a_energies, a_offset, a_productEnergies, 
                a_productMomenta, a_productGains, a_ignoreIncompleteParticles );
}

/* *********************************************************************************************************//**
 * This method modifies the cross section for *this* reaction as
 *
 *          crossSection = a_offset + a_slope * crossSection
 *
 * Either or both of *a_offset* and *a_slope* can be an empty Functions::XYs1d instance or a nullptr.
 * If both *a_offset* and *a_slope* are non-empty Functions::XYs1d instances, the domains of both must be the same.
 * If the returned value is false, no data are changed.
 * 
 * @param       a_offset            [in]    A pointer to a XYs1d function for the offset.
 * @param       a_slope             [in]    A pointer to a XYs1d function for the slope.
 * @param       a_updateMultiGroup  [in]    If true, the multi-group data are also modified.
 *
 * @return                                  true if data can be modified and false otherwise.
 ***********************************************************************************************************/

bool Reaction::modifyCrossSection( Functions::XYs1d const *a_offset, Functions::XYs1d const *a_slope, bool a_updateMultiGroup ) {
/*
    -) FIXME, this does not modify the total cross section for the protare! Does it needed to be?
*/

    ProtareSingle &protare( dynamic_cast<ProtareSingle &>( *root( ) ) );
    Styles::Suite const &styles = protare.styles( );
    Functions::XYs1d const *offset1 = a_offset, *slope1 = a_slope;

    if( a_offset == nullptr ) offset1 = new Functions::XYs1d( );            // Handle nullptr cases.
    if( a_slope == nullptr ) slope1 = new Functions::XYs1d( );

    if( offset1->size( ) == 0 ) {                                           // Handle empty functions.
        if( slope1->size( ) == 0 ) {
            if( a_offset == nullptr ) delete offset1;
            if( a_slope == nullptr ) delete slope1;
            return( false );
        }
        Functions::XYs1d const *offset2 = Functions::XYs1d::makeConstantXYs1d( offset1->axes( ), slope1->domainMin( ), slope1->domainMax( ), 0.0 );
        if( a_offset == nullptr ) delete offset1;
        offset1 = offset2; }
    else if( slope1->size( )  == 0 ) {
        Functions::XYs1d const *slope2 = Functions::XYs1d::makeConstantXYs1d( slope1->axes( ), offset1->domainMin( ), offset1->domainMax( ), 1.0 );
        if( a_slope == nullptr ) delete slope1;
        slope1 = slope2;
    }

    double domainMin = offset1->domainMin( ), domainMax = offset1->domainMax( );

    if( domainMin != slope1->domainMin( ) ) throw Exception( "GIDI::Reaction::modifyCrossSection: offset and slope domainMins differ." );
    if( domainMax != slope1->domainMax( ) ) throw Exception( "GIDI::Reaction::modifyCrossSection: offset and slope domainMaxs differ." );

    Styles::TemperatureInfos temperatureInfos = protare.temperatures( );

    Functions::XYs1d *xys1d = static_cast<Functions::XYs1d *>(
            m_crossSection.findInstanceOfTypeInLineage( styles, temperatureInfos[0].heatedCrossSection( ), GIDI_XYs1dChars ) );
    if( xys1d == nullptr ) throw Exception( "GIDI::Reaction::modifyCrossSection: could not find XYs1d cross section." );

    if( ( xys1d->domainMin( ) >= domainMax ) || ( xys1d->domainMax( ) <= domainMin ) ) {
        if( a_offset == nullptr ) delete offset1;
        if( a_slope == nullptr ) delete slope1;
        return( false );
    }

    double crossSectionDomainMax = xys1d->domainMax( );
    if( xys1d->domainMin( ) > domainMin ) domainMin = xys1d->domainMin( );
    if( crossSectionDomainMax < domainMax ) domainMax = crossSectionDomainMax;

    Functions::XYs1d offset = offset1->domainSlice( domainMin, domainMax, true );
    Functions::XYs1d slope  = slope1->domainSlice( domainMin, domainMax, true );
    if( a_offset == nullptr ) delete offset1;
    if( a_slope == nullptr ) delete slope1;

    for( auto temperatureInfo = temperatureInfos.begin( ); temperatureInfo != temperatureInfos.end( ); ++temperatureInfo ) {
        xys1d = static_cast<Functions::XYs1d *>( 
                m_crossSection.findInstanceOfTypeInLineage( styles, temperatureInfo->heatedCrossSection( ), GIDI_XYs1dChars ) );
        Functions::XYs1d xys1dSliced = xys1d->domainSlice( domainMin, domainMax, true );
        if( xys1dSliced.interpolation( ) != ptwXY_interpolationLinLin ) {
            Functions::XYs1d *temp = xys1dSliced.asXYs1d( true, 1e-3, 1e-6, 1e-6 );
            xys1dSliced = (*temp);
            delete temp;
        }
        Functions::XYs1d modified = offset + slope * xys1dSliced;

        int64_t index1;
        ptwXYPoint *p1;
        ptwXYPoints *ptwXY = xys1d->ptwXY( );
        int64_t length = ptwXY_length( nullptr, ptwXY );
        for( index1 = 0, p1 = ptwXY->points; index1 < length; ++index1, ++p1 ) {
            if( p1->x  < domainMin ) continue;
            if( p1->x >= domainMax ) {
                if( ( p1->x != domainMax ) || ( p1->x != crossSectionDomainMax ) ) break;
            }
            p1->y = modified.evaluate( p1->x );
        }

        Functions::Ys1d *ys1d = m_crossSection.get<Functions::Ys1d>( temperatureInfo->griddedCrossSection( ) );
        std::vector<double> &Ys = ys1d->Ys( );
        Styles::GriddedCrossSection const griddedCrossSection = *styles.get<Styles::GriddedCrossSection>( temperatureInfo->griddedCrossSection( ) );
        nf_Buffer<double> const &grid = griddedCrossSection.grid( ).values( );
        std::size_t start = ys1d->start( ), size = ys1d->size( ), index2;
        for( index2 = 0; index2 < size; ++index2, ++start ) {
            double xValue = grid[start];

            if( xValue  < domainMin ) continue;
            if( xValue >= domainMax ) {
                if( ( xValue != domainMax ) || ( xValue != crossSectionDomainMax ) ) break;
            }
            Ys[index2] = modified.evaluate( xValue );
        }
        if( a_updateMultiGroup ) recalculateMultiGroupData( &protare, *temperatureInfo );
    }

    return( true );
}

/* *********************************************************************************************************//**
 * Thid methid is deprecated, use modifyCrossSection instead. See modifyCrossSection for useage.
 ***********************************************************************************************************/

bool Reaction::modifiedCrossSection( Functions::XYs1d const *a_offset, Functions::XYs1d const *a_slope ) {

    return modifyCrossSection( a_offset, a_slope, false );
}

/* *********************************************************************************************************//**
 * This function recalculates the multi-group data for data labelled with a_temperatureInfo.heatedMultiGroup().
 *
 * @param       a_temperatureInfo   [in]    The temperature for which multi-group data are to be recalculated.
 ***********************************************************************************************************/

void Reaction::recalculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo ) {

    std::string heatedMultiGroupLabel = a_temperatureInfo.heatedMultiGroup( );

    GUPI::Ancestry const *ancestorRoot = root( );
    if( ancestorRoot->moniker( ) != GIDI_topLevelChars ) throw Exception( "Reaction::recalculateMultiGroupData: could not find parent protare." );

    Styles::Suite const &styles = a_protare->styles( );
    if( styles.find( heatedMultiGroupLabel ) == styles.end( ) ) return;

    Styles::HeatedMultiGroup const &heatedMultiGroupStyle = *styles.get<Styles::HeatedMultiGroup const>( heatedMultiGroupLabel );

    std::vector<double> groupBoundaries = heatedMultiGroupStyle.groupBoundaries( a_protare->projectile( ).ID( ) );
    Transporting::MultiGroup multiGroup = Transporting::MultiGroup( "recal", groupBoundaries );

    Transporting::Flux flux( "recal", a_temperatureInfo.temperature( ).value( ) );
    double energies[2] = { groupBoundaries[0], groupBoundaries.back( ) };
    double fluxValues[2] = { 1.0, 1.0 };
    Transporting::Flux_order flux_order( 0, 2, energies, fluxValues );
    flux.addFluxOrder( flux_order );

    MultiGroupCalulationInformation multiGroupCalulationInformation( multiGroup, flux );
    calculateMultiGroupData( a_protare, a_temperatureInfo, heatedMultiGroupLabel, multiGroupCalulationInformation );
}

/* *********************************************************************************************************//**
 * This methods calculates multi-group data for all needed components and adds each component's multi-group with label *a_heatedMultiGroupLabel*.
 *
 * @param   a_temperatureInfo                   [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param   a_heatedMultiGroupLabel             [in]    The label of the style for the multi-group data being added.
 * @param   a_multiGroupCalulationInformation   [in]    Store multi-group boundary and flux data used for multi-grouping.
 * @param   a_crossSectionXYs1d                 [in[    The cross section weight.
 ***********************************************************************************************************/

// FIXME maybe, as upscatter is currently not handled.
    
void Reaction::calculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo, 
                std::string const &a_heatedMultiGroupLabel, MultiGroupCalulationInformation const &a_multiGroupCalulationInformation ) {

    Styles::Suite const &styles = a_protare->styles( );
    Transporting::MultiGroup multiGroup = a_multiGroupCalulationInformation.m_multiGroup;

    std::vector<double> flatValues { multiGroup[0], 1.0, multiGroup[multiGroup.size( ) - 1], 1.0 };
    Functions::XYs1d flatFunction( Axes( ), ptwXY_interpolationLinLin, flatValues );

    Functions::XYs1d const *crossSectionXYs1d = static_cast<Functions::XYs1d *>( 
                m_crossSection.findInstanceOfTypeInLineage( styles, a_temperatureInfo.heatedCrossSection( ), GIDI_XYs1dChars ) );

    calculate1dMultiGroupDataInComponent( a_protare, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, m_crossSection, flatFunction );
    calculate1dMultiGroupDataInComponent( a_protare, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, m_availableEnergy, *crossSectionXYs1d );
    calculate1dMultiGroupDataInComponent( a_protare, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, m_availableMomentum, *crossSectionXYs1d );

    if( m_outputChannel != nullptr ) 
        m_outputChannel->calculateMultiGroupData( a_protare, a_temperatureInfo, a_heatedMultiGroupLabel, a_multiGroupCalulationInformation, *crossSectionXYs1d );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Reaction::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {
    
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string attributes;

    attributes += a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
    attributes += a_writeInfo.addAttribute( GIDI_ENDF_MT_Chars, intToString( ENDF_MT( ) ) );
    if( m_fissionGenre != "" ) attributes += a_writeInfo.addAttribute( GIDI_fissionGenreChars, m_fissionGenre );
    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    m_doubleDifferentialCrossSection.toXMLList( a_writeInfo, indent2 );    
    m_crossSection.toXMLList( a_writeInfo, indent2 );    
    if( m_outputChannel != nullptr ) m_outputChannel->toXMLList( a_writeInfo, indent2 );
    m_availableEnergy.toXMLList( a_writeInfo, indent2 );
    m_availableMomentum.toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeEnder( moniker( ) );
}

/* *********************************************************************************************************//**
 * Calculates the ENDL C and S values for a ENDF MT value.
 *
 * @param ENDF_MT   [in]    The ENDF MT value.
 * @param ENDL_C    [out]   The ENDL C value for the ENDF MT value.
 * @param ENDL_S    [out]   The ENDL S value for the ENDF MT value.
 *
 * @return                  Returns 0 if the ENDF MT value is valid and 1 otherwise.
 ***********************************************************************************************************/

int ENDL_CFromENDF_MT( int ENDF_MT, int *ENDL_C, int *ENDL_S ) {

    int MT1_50ToC[] = { 1,   10,  -3,   11,   -5,    0,    0,    0,    0,  -10,
                       32,    0,   0,    0,    0,   12,   13,   15,   15,   15,
                       15,   26,  36,   33,  -25,    0,  -27,   20,   27,  -30,
                        0,   22,  24,   25,  -35,  -36,   14,   15,    0,    0,
                       29,   16,   0,   17,   34,    0,    0,    0,    0 };
    int MT100_200ToC[] = { -101,   46,   40,   41,   42,   44,   45,   37, -109,    0,
                             18,   48, -113, -114,   19,   39,   47,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                              0, -152, -153, -154,   43, -156, -157,   23,   31, -160,
                           -161, -162, -163, -164, -165, -166, -167, -168, -169, -170,
                           -171, -172, -173, -174, -175, -176, -177, -178, -179, -180,
                           -181, -182, -183, -184, -185, -186, -187, -188,   28, -190,
                           -191, -192,   38, -194, -195, -196, -197, -198, -199, -200 };
    *ENDL_C = 0;
    *ENDL_S = 0;
    if( ENDF_MT <= 0 ) {
        *ENDL_C = -ENDF_MT;
        return( 1 );
    }
    if( ENDF_MT > 1572 ) return( 1 );
    if( ENDF_MT < 50 ) {
        *ENDL_C = MT1_50ToC[ENDF_MT - 1]; }
    else if( ENDF_MT <= 91 ) {
        *ENDL_C = 11;
        if( ENDF_MT != 91 ) *ENDL_S = 1; }
    else if( ( ENDF_MT > 100 ) && ( ENDF_MT <= 200 ) ) {
        *ENDL_C = MT100_200ToC[ENDF_MT - 101]; }
    else if( ( ENDF_MT == 452 ) || ( ENDF_MT == 455 ) || ( ENDF_MT == 456 ) || ( ENDF_MT == 458 ) ) {
        *ENDL_C = 15;
        if( ENDF_MT == 455 ) *ENDL_S = 7; }
    else if( ( ENDF_MT >= 502 ) && ( ENDF_MT <= 572 ) ) {
        if( ENDF_MT == 502 ) {
            *ENDL_C = 71; }
        else if( ENDF_MT == 504 ) {
            *ENDL_C = 72; }
        else if( ( ENDF_MT >= 515 ) && ( ENDF_MT <= 517 ) ) {
            *ENDL_C = 74; }
        else if( ENDF_MT == 522 ) {
            *ENDL_C = 73; }
        else if( ( ENDF_MT >= 534 ) && ( ENDF_MT <= 572 ) ) {
            *ENDL_C = 73; } }
    else if( ENDF_MT >= 600 ) {
        if( ENDF_MT < 650 ) {
            *ENDL_C = 40;
            if( ENDF_MT != 649 ) *ENDL_S = 1; }
        else if( ENDF_MT < 700 ) {
            *ENDL_C = 41;
            if( ENDF_MT != 699 ) *ENDL_S = 1; }
        else if( ENDF_MT < 750 ) {
            *ENDL_C = 42;
            if( ENDF_MT != 749 ) *ENDL_S = 1; }
        else if( ENDF_MT < 800 ) {
            *ENDL_C = 44;
            if( ENDF_MT != 799 ) *ENDL_S = 1; }
        else if( ENDF_MT < 850 ) {
            *ENDL_C = 45;
            if( ENDF_MT != 849 ) *ENDL_S = 1; }
        else if( ( ENDF_MT >= 875 ) && ( ENDF_MT <= 891 ) ) {
            *ENDL_C = 12;
            if( ENDF_MT != 891 ) *ENDL_S = 1; }
        else if( ( ENDF_MT >= 1534 ) && (ENDF_MT <= 1572 ) ) {
            *ENDL_C = 72;
        }
    }

    return( 0 );
}

}
