/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <algorithm>

#include "GIDI.hpp"

namespace GIDI {

/*! \class ProtareComposite
 * Class to store a list of GNDS <**reactionSuite**> node.
 */

/* *********************************************************************************************************//**
 * ProtareComposite constructor.
 *
 * @param a_construction        [in]     Used to pass user options to the constructor.
 ***********************************************************************************************************/

ProtareComposite::ProtareComposite( LUPI_maybeUnused Construction::Settings const &a_construction ) {

}

/* *********************************************************************************************************//**
 ******************************************************************/

ProtareComposite::~ProtareComposite( ) {

    for( std::vector<Protare *>::const_iterator iter = m_protares.begin( ); iter < m_protares.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Appends *a_protare* to the list of Protares.
 *
 * @param a_protare     [in]    The Protare to add to *this* instance.
 ******************************************************************/

void ProtareComposite::append( Protare *a_protare ) {

    m_protares.push_back( a_protare );
}

/* *********************************************************************************************************//**
 * Returns the GNDS format version for the (a_index+1)^th Protare.
 *
 * @param           a_index [in]    The index of the Protare whose format version is returned.
 *
 * @return                          The format version.
 ******************************************************************/

LUPI::FormatVersion const &ProtareComposite::formatVersion( std::size_t a_index ) const {

    return( m_protares[a_index]->formatVersion( ) );
}

/* *********************************************************************************************************//**
 * Returns the file name for the (a_index+1)^th Protare.
 *
 * @param           a_index [in]    The index of the Protare whose file name is returned.
 *
 * @return                          The file name.
 ******************************************************************/

std::string const &ProtareComposite::fileName( std::size_t a_index ) const {
    
    return( m_protares[a_index]->fileName( ) );
}

/* *********************************************************************************************************//**
 * Returns the real file name for the (a_index+1)^th Protare.
 *
 * @param           a_index [in]    The index of the Protare whose real file name is returned.
 *
 * @return                          The real file name.
 ******************************************************************/

std::string const &ProtareComposite::realFileName( std::size_t a_index ) const {
    
    return( m_protares[a_index]->realFileName( ) );
}

/* *********************************************************************************************************//**
 * Returns the list of libraries for the (a_index+1)^th contained Protare.
 * 
 * @param           a_index     [in]    The index of the Protare whose libraries are to be returned.
 *
 * @return                              The list of libraries.
 ******************************************************************/

std::vector<std::string> ProtareComposite::libraries( std::size_t a_index ) const {

    return( m_protares[a_index]->libraries( ) );
}

/* *********************************************************************************************************//**
 * Returns the evaluation for the (a_index+1)^th Protare.
 *
 * @param           a_index [in]    The index of the Protare whose evaluation is returned.
 *
 * @return                          The evaluation.
 ******************************************************************/

std::string const &ProtareComposite::evaluation( std::size_t a_index ) const {
    
    return( m_protares[a_index]->evaluation( ) );
}

/* *********************************************************************************************************//**
 * Returns the projectile frame for the (a_index+1)^th Protare.
 *
 * @param           a_index [in]    The index of the Protare whose projectile frame is returned.
 *
 * @return                          The projectile frame.
 ******************************************************************/

Frame ProtareComposite::projectileFrame( std::size_t a_index ) const {
    
    return( m_protares[a_index]->projectileFrame( ) );
}

/* *********************************************************************************************************//**
 * Returns the number of **ProtareSingle**s contained in *this*.
 *
 * @return                              Returns the number of contained **ProtareSingle**s..
 ***********************************************************************************************************/

std::size_t ProtareComposite::numberOfProtares( ) const {

    std::size_t number = 0;

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) number += m_protares[i1]->numberOfProtares( );

    return( number );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the (a_index - 1)th **ProtareSingle**.
 *
 * @param a_index               [in]    Index of the **ProtareSingle** to return.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

ProtareSingle *ProtareComposite::protare( std::size_t a_index ) {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        std::size_t number = m_protares[i1]->numberOfProtares( );

        if( number > a_index ) return( m_protares[i1]->protare( a_index ) );
        a_index -= number;
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the (a_index - 1)th **ProtareSingle**.
 *
 * @param a_index               [in]    Index of the **ProtareSingle** to return.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

ProtareSingle const *ProtareComposite::protare( std::size_t a_index ) const {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        std::size_t number = m_protares[i1]->numberOfProtares( );

        if( number > a_index ) return( m_protares[i1]->protare( a_index ) );
        a_index -= number;
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the number of LazyParsingHelperForms instantiated.
 *
 * @return              The number of LazyParsingHelperForms instantiated.
 ******************************************************************/

int ProtareComposite::numberOfLazyParsingHelperForms( ) const {

    int numberOfLazyParsingHelperForms1 = 0;

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) numberOfLazyParsingHelperForms1 += m_protares[i1]->numberOfLazyParsingHelperForms( );

    return( numberOfLazyParsingHelperForms1 );
}

/* *********************************************************************************************************//**
 * Returns the number of instantiated LazyParsingHelperForms replaced with the appropriate form.
 *
 * @return              The number of LazyParsingHelperForms replaced.
 ******************************************************************/

int ProtareComposite::numberOfLazyParsingHelperFormsReplaced( ) const {

    int numberOfLazyParsingHelperFormsReplaced1 = 0;

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) numberOfLazyParsingHelperFormsReplaced1 += m_protares[i1]->numberOfLazyParsingHelperFormsReplaced( );

    return( numberOfLazyParsingHelperFormsReplaced1 );
}

/* *********************************************************************************************************//**
 * Returns a threshold factor for the projectile hitting the target.
 *
 * @return              The threshold factor.
 ******************************************************************/

double ProtareComposite::thresholdFactor( ) const {
    
    return( m_protares[0]->thresholdFactor( ) );
}

/* *********************************************************************************************************//**
 * Returns the Documentation_1_10::Suite from the first protare in *m_protares*.
 *
 * @return              The Documentation_1_10::Suite.
 ******************************************************************/

Documentation_1_10::Suite &ProtareComposite::documentations( ) {

    return( m_protares[0]->documentations( ) );
}

/* *********************************************************************************************************//**
 * Returns the style with label *a_label* from the first Protare in *m_protares*.
 *
 * @param           a_label     [in]    The label of the requested style.
 * @return                              The style with label *a_label*.
 ******************************************************************/

Styles::Base &ProtareComposite::style( std::string const &a_label ) {

    return( m_protares[0]->style( a_label ) );
}

/* *********************************************************************************************************//**
 * Returns the Styles::Suite from the first Protare in *m_protares*.
 *  
 * @return              The Styles::Suite.
 ******************************************************************/

Styles::Suite &ProtareComposite::styles( ) {

    return( m_protares[0]->styles( ) );
}

/* *********************************************************************************************************//**
 * Returns the Styles::Suite from the first Protare in *m_protares*.
 *
 * @return              The Styles::Suite.
 ******************************************************************/

Styles::Suite const &ProtareComposite::styles( ) const {

    return( m_protares[0]->styles( ) );
}

/* *********************************************************************************************************//**
 * Returns the intid for the requested particle or -1 if the particle is not in *m_protare* PoPs database.
 *
 * @param a_id                 [in]    The GNDS PoPs id for particle whose intd is requested.
 *
 * @return                             C++ int for the requested particle or -1 if particle is not in PoPs. 
 ******************************************************************/

int ProtareComposite::intid( std::string const &a_id ) const {

    int intid1 = -1;

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        intid1 = m_protares[i1]->intid( a_id );
        if( intid1 > -1 ) break;
    }

    return( intid1 );
}

/* *********************************************************************************************************//**
 * Calls productIDs for each Protare contained in *this*.
 *
 * @param a_ids                 [in]    The unique list of product indices.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_transportablesOnly  [in]    If *true* only transportable particle ids are added to *a_ids*.
 ******************************************************************/

void ProtareComposite::productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) m_protares[i1]->productIDs( a_ids, a_particles, a_transportablesOnly );
}

/* *********************************************************************************************************//**
 * Determines the maximum Legredre order present in the multi-group transfer matrix for a give product for a give label.
 * Loops over all contained Protares to determine the maximum Legredre order.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    The id of the requested product.
 *
 * @return                          The maximum Legredre order. If no transfer matrix data are present for the requested product, -1 is returned.
 ***********************************************************************************************************/

int ProtareComposite::maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                    Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    int maximumLegendreOrder1 = -1;

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        int maximumLegendreOrder2 = m_protares[i1]->maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID );

        if( maximumLegendreOrder2 > maximumLegendreOrder1 ) maximumLegendreOrder1 = maximumLegendreOrder2;
    }

    return( maximumLegendreOrder1 );
}

/* *********************************************************************************************************//**
 * Returns a list of all process temperature data. For each temeprature, the labels for its
 *
 *   + heated cross section data,
 *   + gridded cross section data,
 *   + multi-group data, and
 *   + multi-group upscatter data.
 *
 * are returned. If no data are present for a give data type (e.g., gridded cross section, multi-group upscatter), its label is an empty std::string.
 *
 * @return  The list of temperatures and their labels via an Styles::TemperatureInfos instance. The Styles::TemperatureInfos class
 *          has many (if not all) the method of a std::vector.
 ***********************************************************************************************************/

Styles::TemperatureInfos ProtareComposite::temperatures( ) const {

    return( m_protares[0]->temperatures( ) );
}

/* *********************************************************************************************************//**
 * Returns the number of reactions for all Protares contained in *this*.
 *  
 * @return              The total number of reactions.
 ******************************************************************/

std::size_t ProtareComposite::numberOfReactions( ) const {

    std::size_t numberOfReactions1 = 0;

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) numberOfReactions1 += m_protares[i1]->numberOfReactions( );

    return( numberOfReactions1 );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th reaction.
 *
 * @param a_index           [in]    The index of the requested reaction.
 * @return                          The (*a_index*+1)th reaction.
 ***********************************************************************************************************/

Reaction *ProtareComposite::reaction( std::size_t a_index ) {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        std::size_t numberOfReactions1 = m_protares[i1]->numberOfReactions( );

        if( a_index < numberOfReactions1 ) return( m_protares[i1]->reaction( a_index ) );
        a_index -= numberOfReactions1;
    }

    throw Exception( "ProtareComposite::reaction: index out of range" );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th reaction.
 *
 * @param a_index           [in]    The index of the requested reaction.
 * @return                          The (*a_index*+1)th reaction.
 ***********************************************************************************************************/

Reaction const *ProtareComposite::reaction( std::size_t a_index ) const {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        std::size_t numberOfReactions1 = m_protares[i1]->numberOfReactions( );

        if( a_index < numberOfReactions1 ) return( m_protares[i1]->reaction( a_index ) );
        a_index -= numberOfReactions1;
    }

    throw Exception( "ProtareComposite::reaction: index out of range" );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th reaction or **nullptr**. If the indexed reaction is deactivated or exlucded, 
 * a **nullptr** is returned.
 *
 * @param a_index               [in]    The index of the requested reaction.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the cross section.
 *
 * @return                          The (*a_index*+1)th reaction or **nullptr**.
 ***********************************************************************************************************/

Reaction const *ProtareComposite::reaction( std::size_t a_index, Transporting::MG const &a_settings, 
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        std::size_t numberOfReactions1 = m_protares[i1]->numberOfReactions( );

        if( a_index < numberOfReactions1 ) return( m_protares[i1]->reaction( a_index, a_settings, a_reactionsToExclude ) );
        a_index -= numberOfReactions1;
    }

    throw Exception( "ProtareComposite::reaction: index out of range" );
}

/* *********************************************************************************************************//**
 * Returns the number of orphanProduct for all Protares contained in *this*.
 *  
 * @return              The total number of orphanProducts.
 ******************************************************************/
    
std::size_t ProtareComposite::numberOfOrphanProducts( ) const {

    std::size_t numberOfOrphanProducts1 = 0;

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) numberOfOrphanProducts1 += m_protares[i1]->numberOfOrphanProducts( );

    return( numberOfOrphanProducts1 );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th orphanProduct.
 *
 * @param a_index           [in]    The index of the requested orphanProduct.
 * @return                          The (*a_index*+1)th orphanProduct.
 ***********************************************************************************************************/

Reaction *ProtareComposite::orphanProduct( std::size_t a_index ) {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        std::size_t numberOfOrphanProducts1 = m_protares[i1]->numberOfOrphanProducts( );
        
        if( a_index < numberOfOrphanProducts1 ) return( m_protares[i1]->orphanProduct( numberOfOrphanProducts1 ) );
        a_index -= numberOfOrphanProducts1;
    }

    throw Exception( "ProtareComposite::orphanProduct: index out of range" );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th orphanProduct.
 *
 * @param a_index           [in]    The index of the requested orphanProduct.
 * @return                          The (*a_index*+1)th orphanProduct.
 ***********************************************************************************************************/

Reaction const *ProtareComposite::orphanProduct( std::size_t a_index ) const {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        std::size_t numberOfOrphanProducts1 = m_protares[i1]->numberOfOrphanProducts( );

        if( a_index < numberOfOrphanProducts1 ) return( m_protares[i1]->orphanProduct( numberOfOrphanProducts1 ) );
        a_index -= numberOfOrphanProducts1;
    }

    throw Exception( "ProtareComposite::orphanProduct: index out of range" );
}

/* *********************************************************************************************************//**
 * Re-indexs the reactions in the reactions, orphanProducts and fissionComponents suites.
 *
 ***********************************************************************************************************/
    
void ProtareComposite::updateReactionIndices( LUPI_maybeUnused std::size_t a_offset ) const {

    std::size_t reactionOffset = 0;

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        m_protares[i1]->updateReactionIndices( reactionOffset );
        reactionOffset += m_protares[i1]->numberOfReactions( );
    }
}

/* *********************************************************************************************************//**
 * Returns true if at least one reaction contains a fission channel.
 *
 * @return      true if at least one reaction contains a fission channel and false otherwise.
 ***********************************************************************************************************/

bool ProtareComposite::hasFission( ) const {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {

        if( m_protares[i1]->hasFission( ) ) return( true );
    }
    return( false );
}

/* *********************************************************************************************************//**
 * Returns **false* if protare has delayed fission neutrons for an active reaction and they are not complete; otherwise, returns **true**.
 *
 * @return      bool
 ***********************************************************************************************************/

bool ProtareComposite::isDelayedFissionNeutronComplete( ) const {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
    
        if( !m_protares[i1]->isDelayedFissionNeutronComplete( ) ) return( true );
    }

    return( true );
}

/* *********************************************************************************************************//**
 * Returns the multi-group boundaries for the requested label and product.
 *
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID       [in]    ID for the requested product.
 *
 * @return                          List of multi-group boundaries.
 ***********************************************************************************************************/

std::vector<double> ProtareComposite::groupBoundaries( Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    return( m_protares[0]->groupBoundaries( a_settings, a_temperatureInfo, a_productID ) );
}

/* *********************************************************************************************************//**
 * Returns the inverse speeds for the requested label. The label must be for a heated multi-group style.
 *
 * @param a_smr             [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings        [in]    Specifies the requested label.
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 *
 * @return                          List of inverse speeds.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupInverseSpeed( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo ) const {

    return( m_protares[0]->multiGroupInverseSpeed( a_smr, a_settings, a_temperatureInfo ) );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total cross section for the requested label. This is summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the cross section.
 * @param a_label               [in]    If not an empty string, this is used as the label for the form to return and the *a_temperatureInfo* labels are ignored.
 *
 * @return                              The requested multi-group cross section as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupCrossSection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude, std::string const &a_label ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupCrossSection( a_smr, a_settings, a_temperatureInfo, a_reactionsToExclude, a_label );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total Q for the requested label. This is a cross section weighted multiplicity
 * summed over all reactions
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_final               [in]    If false, only the Q for the primary reactions are return, otherwise, the Q for the final reactions.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the Q.
 *
 * @return                              The requested multi-group Q as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, bool a_final, bool a_effectivePhotoAtomic,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final, a_effectivePhotoAtomic, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total multiplicity for the requested label for the requested product. This is a cross section weighted multiplicity.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    Id for the requested product.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the multiplicity.
 *
 * @return                              The requested multi-group multiplicity as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total fission neutron multiplicity for the requested label. This is a cross section weighted multiplicity.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the multiplicity.
 *
 * @return                              The requested multi-group fission neutron multiplicity as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupFissionNeutronMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupFissionNeutronMultiplicity( a_smr, a_settings, a_temperatureInfo, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total fission gamma multiplicity for the requested label. This is a cross section weighted multiplicity.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the multiplicity.
 *
 * @return                              The requested multi-group fission neutron multiplicity as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupFissionGammaMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupFissionGammaMultiplicity( a_smr, a_settings, a_temperatureInfo, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total product matrix for the requested label for the requested product id for the requested Legendre order.
 * If no data are found, an empty GIDI::Matrix is returned.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label and if delayed neutrons should be included.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_productID           [in]    PoPs id for the requested product.
 * @param a_order               [in]    Requested product matrix, Legendre order.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the product matrix.
 *
 * @return                              The requested multi-group product matrix as a GIDI::Matrix.
 ***********************************************************************************************************/

Matrix ProtareComposite::multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, 
                std::string const &a_productID, std::size_t a_order, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Matrix matrix( 0, 0 );
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        matrix += m_protares[i1]->multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order,
                excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( matrix );
}

/* *********************************************************************************************************//**
 * Like ProtareComposite::multiGroupProductMatrix, but only returns the fission neutron, transfer matrix.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label and if delayed neutrons should be included.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_order               [in]    Requested product matrix, Legendre order.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the fission matrix.
 *
 * @return                              The requested multi-group neutron fission matrix as a GIDI::Matrix.
 ***********************************************************************************************************/

Matrix ProtareComposite::multiGroupFissionMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::size_t a_order,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Matrix matrix( 0, 0 );
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        matrix += m_protares[i1]->multiGroupFissionMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_order, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( matrix );
}

/* *********************************************************************************************************//**
 * Returns the multi-group transport correction for the requested label. The transport correction is calculated from the transfer matrix
 * for the projectile id for the Legendre order of *a_order + 1*.
 *
 * @param a_smr                     [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings                [in]    Specifies the requested label.
 * @param a_temperatureInfo         [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles               [in]    The list of particles to be transported.
 * @param a_order                   [in]    Maximum Legendre order for transport. The returned transport correction is for the next higher Legender order.
 * @param a_transportCorrectionType [in]    Requested transport correction type.
 * @param a_temperature             [in]    The temperature of the flux to use when collapsing. Pass to the GIDI::collapse method.
 * @param a_reactionsToExclude      [in]    A list of reaction indices that are to be ignored when calculating the transport correction.
 *
 * @return                                  The requested multi-group transport correction as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupTransportCorrection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::size_t a_order, 
                TransportCorrectionType a_transportCorrectionType, double a_temperature, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupTransportCorrection( a_smr, a_settings, a_temperatureInfo, a_particles, a_order, 
                a_transportCorrectionType, a_temperature, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total available energy for the requested label. This is a cross section weighted available energy
 * summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the available energy.
 *
 * @return                              The requested multi-group available energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupAvailableEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 )  {
        vector += m_protares[i1]->multiGroupAvailableEnergy( a_smr, a_settings, a_temperatureInfo, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total average energy for the requested label for the requested product. This is a cross section weighted average energy
 * summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    Particle id for the requested product.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the average energy.
 *
 * @return                              The requested multi-group average energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total deposition energy for the requested label. This is a cross section weighted deposition energy
 * summed over all reactions. The deposition energy is calculated by subtracting the average energy from each transportable particle
 * from the available energy. The list of transportable particles is specified via the list of particle specified in the *a_settings* argument.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label and the products that are transported.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the deposition energy.
 *
 * @return                              The requested multi-group deposition energy as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupDepositionEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupDepositionEnergy( a_smr, a_settings, a_temperatureInfo, a_particles, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total available momentum for the requested label. This is a cross section weighted available momentum
 * summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the available momentum.
 *
 * @return                              The requested multi-group available momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupAvailableMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 )  {
        vector += m_protares[i1]->multiGroupAvailableMomentum( a_smr, a_settings, a_temperatureInfo, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total average momentum for the requested label for the requested product. This is a cross section weighted average momentum
 * summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    Particle id for the requested product.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the average momentum.
 *
 * @return                              The requested multi-group average momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, total deposition momentum for the requested label. This is a cross section weighted deposition momentum
 * summed over all reactions. The deposition momentum is calculated by subtracting the average momentum from each transportable particle
 * from the available momentum. The list of transportable particles is specified via the list of particle specified in the *a_settings* argument.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_particles           [in]    The list of particles to be transported.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the deposition momentum.
 *
 * @return                              The requested multi-group deposition momentum as a GIDI::Vector.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupDepositionMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupDepositionMomentum( a_smr, a_settings, a_temperatureInfo, a_particles, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 * Returns the multi-group, gain for the requested particle and label. This is a cross section weighted gain summed over all reactions.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    The PoPs' id for the particle whose gain is to be calculated.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the gain.
 *
 * @return                              The requested multi-group gain as a **GIDI::Vector**.
 ***********************************************************************************************************/

Vector ProtareComposite::multiGroupGain( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    Vector vector;
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        vector += m_protares[i1]->multiGroupGain( a_smr, a_settings, a_temperatureInfo, a_productID, excludeReactionsSet );
        excludeReactionsSetAdjust( excludeReactionsSet, *m_protares[i1] );
    }

    return( vector );
}

/* *********************************************************************************************************//**
 *
 *
 * @return      A list of label, mu cutoff pairs.
 ***********************************************************************************************************/

stringAndDoublePairs ProtareComposite::muCutoffForCoulombPlusNuclearElastic( ) const {

    stringAndDoublePairs muCutoffs;

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        stringAndDoublePairs muCutoffs2 = m_protares[i1]->muCutoffForCoulombPlusNuclearElastic( );

        for( stringAndDoublePairs::iterator iter = muCutoffs2.begin( ); iter < muCutoffs2.end( ); ++iter ) muCutoffs.push_back( *iter );
    }

    return( muCutoffs );
}

/* *********************************************************************************************************//**
 * Returns the list of DelayedNeutronProduct instances.
 * 
 * @return                      The list of delayed neutrons.
 ***********************************************************************************************************/
 
DelayedNeutronProducts ProtareComposite::delayedNeutronProducts( ) const {

    DelayedNeutronProducts delayedNeutronProducts1;

    if( hasFission( ) ) {
        for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
            DelayedNeutronProducts delayedNeutronProducts2 = m_protares[i1]->delayedNeutronProducts( );

            delayedNeutronProducts1.insert( delayedNeutronProducts1.end( ), delayedNeutronProducts2.begin( ), delayedNeutronProducts2.end( ) );
        }
    }

    return( delayedNeutronProducts1 );
}

/* *********************************************************************************************************//**
 * Calls the **incompleteParticles** method for each **ProtareSingle** in *this*.
 *
 * @param a_settings        [in]    Specifies the requested label.
 * @param       a_incompleteParticles   [out]   The list of particles whose **completeParticle** method returns *false*.
 ***********************************************************************************************************/

void ProtareComposite::incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const {

    for( std::size_t i1 = 0; i1 < m_protares.size( ); ++i1 ) {
        m_protares[i1]->incompleteParticles( a_settings, a_incompleteParticles );
    }
}

}
