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

/*! \class ProtareTNSL
 * Class to store <**reactionSuite**> nodes required for thermal neutron scattering law.
 */

/* *********************************************************************************************************//**
 * ProtareTNSL constructor.
 *
 * @param a_construction        [in]     Used to pass user options to the constructor.
 * @param a_protare             [in]     The non-TNSL **ProtareSingle**.
 * @param a_TNSL                [in]     The TNSL **ProtareSingle**.
 ***********************************************************************************************************/

ProtareTNSL::ProtareTNSL( LUPI_maybeUnused Construction::Settings const &a_construction, ProtareSingle *a_protare, ProtareSingle *a_TNSL ) :
        m_protare( a_protare ),
        m_TNSL( a_TNSL ),
        m_elasticReaction( nullptr ) {

    if( a_protare->projectile( ).ID( ) != PoPI::IDs::neutron ) throw Exception( "ProtareTNSL::ProtareTNSL: a_protare neutron as target." );
    if( a_TNSL->projectile( ).ID( ) != PoPI::IDs::neutron ) throw Exception( "ProtareTNSL::ProtareTNSL: a_TNSL not thermal neutron scattering protare." );

    setProjectile( a_protare->projectile( ) );
    setTarget( a_protare->target( ) );

    for( std::size_t i1 = 0; i1 < m_protare->numberOfReactions( ); ++i1 ) {
        Reaction *reaction = m_protare->reaction( i1 );

        if( reaction->ENDF_MT( ) == 2 ) {
            m_elasticReaction = reaction;
            break;
        }
    }

    if( m_elasticReaction == nullptr ) throw Exception( "ProtareTNSL::ProtareTNSL: could not find elastic reaction in a_protare." );

    m_maximumTNSL_MultiGroupIndex[""] = 0;                             // In case temperatureInfo.heatedMultiGroup( ) is an empty string. Is this needed?
    Styles::TemperatureInfos temperatures = m_TNSL->temperatures( );
    Transporting::MG settings( m_TNSL->projectile( ).ID( ), Transporting::Mode::multiGroup, Transporting::DelayedNeutrons::off );
    LUPI::StatusMessageReporting smr1;
    for( auto iter = temperatures.begin( ); iter != temperatures.end( ); ++iter ) {
        std::string label = iter->heatedMultiGroup( );

        if( label != "" ) {
            std::size_t i1 = 0;
            std::size_t maximumMultiGroupIndex = 0;

            Vector crossSection = m_TNSL->multiGroupCrossSection( smr1, settings, *iter );
            smr1.clear( );
            for( ; i1 < crossSection.size( ); ++i1 ) {
                if( crossSection[i1] == 0 ) break;
            }
            maximumMultiGroupIndex = i1;
                                                                    // The check of multiplicity is needed because FUDGE has a inconsistency in what maximum energy 
                                                                    // to use for TNSL calculations. The following lines can be remove when this issue is resolved in FUDGE.
                                                                    // This should not be needed for data produced after 15/Jan/2021 or probably earlier as FUDGE was fixed.
            Vector multiplicity = m_TNSL->multiGroupMultiplicity( smr1, settings, *iter, PoPI::IDs::neutron );
            smr1.clear( );
            for( i1 = 0; i1 < multiplicity.size( ); ++i1 ) {
                if( multiplicity[i1] == 0 ) break;
            }
            if( i1 < maximumMultiGroupIndex ) maximumMultiGroupIndex = i1;

            m_maximumTNSL_MultiGroupIndex[label] = maximumMultiGroupIndex;
        }
    }

    m_elasticReaction->modifiedMultiGroupElasticForTNSL( m_maximumTNSL_MultiGroupIndex );

// FIXME do I need to check that data are consistence.
}

/* *********************************************************************************************************//**
 ******************************************************************/

ProtareTNSL::~ProtareTNSL( ) {

    delete m_protare;
    delete m_TNSL;
}

/* *********************************************************************************************************//**
 * Returns the maximum number of usable multi-groups for the thermal neutron scattering law protare for the request multi-group label in *a_temperatureInfo*.
 *
 * @param a_temperatureInfo [in]    Specifies the temperature and labels use to lookup the requested data.
 *
 * @return                          The last multigroup index for which the TNSL protare has cross section data. Above this index, the cross section, et al. data must come from the standard protare.
 ******************************************************************/

std::size_t ProtareTNSL::maximumTNSL_MultiGroupIndex( Styles::TemperatureInfo const &a_temperatureInfo ) const {

    return( m_maximumTNSL_MultiGroupIndex.at( a_temperatureInfo.heatedMultiGroup( ) ) );
}

/* *********************************************************************************************************//**
 * Removes the elastic component from *a_vector* and adds in the TNSL component.
 *
 * @param a_settings            [in]        Specifies the requested label.
 * @param a_temperatureInfo     [in]        Specifies the temperature and labels use to lookup the requested data.
 * @param a_vector              [in/out]    The vector from the non TNSL protare.
 * @param a_vectorElastic       [in]        The vector from the elastic reactino from the non TNSL protare.
 * @param a_vectorTNSL          [in]        The vector from the TNSL protare.
 ******************************************************************/

void ProtareTNSL::combineVectors( LUPI_maybeUnused Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, Vector &a_vector, Vector const &a_vectorElastic, Vector const &a_vectorTNSL ) const {

    if( a_vectorTNSL.size( ) == 0 ) return;

    std::size_t maximumMultiGroupIndex = m_maximumTNSL_MultiGroupIndex.at( a_temperatureInfo.heatedMultiGroup( ) );

    for( std::size_t i1 = 0; i1 < maximumMultiGroupIndex; ++i1 ) a_vector[i1] += a_vectorTNSL[i1] - a_vectorElastic[i1];
}

/* *********************************************************************************************************//**
 * Removes the elastic component from *a_matrix* and adds in the TNSL component.
 *
 * @param a_settings            [in]        Specifies the requested label.
 * @param a_temperatureInfo     [in]        Specifies the temperature and labels use to lookup the requested data.
 * @param a_matrix              [in/out]    The matrix from the non TNSL protare.
 * @param a_matrixElastic       [in]        The matrix from the elastic reactino from the non TNSL protare.
 * @param a_matrixTNSL          [in]        The matrix from the TNSL protare.
 ******************************************************************/

void ProtareTNSL::combineMatrices( LUPI_maybeUnused Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, Matrix &a_matrix, Matrix const &a_matrixElastic, Matrix const &a_matrixTNSL ) const {

    if( a_matrixTNSL.size( ) == 0 ) return;

    std::size_t maximumMultiGroupIndex = m_maximumTNSL_MultiGroupIndex.at( a_temperatureInfo.heatedMultiGroup( ) );

    for( std::size_t i1 = 0; i1 < maximumMultiGroupIndex; ++i1 ) {
        Vector const &rowTNSL = a_matrixTNSL[i1];
        Vector const &rowElastic = a_matrixElastic[i1];
        Vector &row = a_matrix[i1];

        for( std::size_t i2 = 0; i2 < a_matrixTNSL.size( ); ++i2 ) row[i2] += rowTNSL[i2] - rowElastic[i2];
    }
}

/* *********************************************************************************************************//**
 * Returns the GNDS format version for the (a_index+1)^th Protare. The index **a_index** can only be 0 (normal protare) or 1 (TNSL protare).
 *
 * @param           a_index [in]    The index of the Protare whose format version is returned.
 *
 * @return                          The format version.
 ******************************************************************/

LUPI::FormatVersion const &ProtareTNSL::formatVersion( std::size_t a_index ) const {

    if( a_index == 0 ) return( m_protare->formatVersion( ) );
    if( a_index == 1 ) return( m_TNSL->formatVersion( ) );
    throw Exception( "ProtareTNSL::formatVersion: index can only be 0 or 1." );
}

/* *********************************************************************************************************//**
 * Returns the file name for the (a_index+1)^th Protare. The index **a_index** can only be 0 (normal protare) or 1 (TNSL protare).
 *
 * @param           a_index [in]    The index of the Protare whose file name is returned.
 *
 * @return                          The file name.
 ******************************************************************/

std::string const &ProtareTNSL::fileName( std::size_t a_index ) const {

    if( a_index == 0 ) return( m_protare->fileName( ) );
    if( a_index == 1 ) return( m_TNSL->fileName( ) );
    throw Exception( "ProtareTNSL::fileName: index can only be 0 or 1." );
}

/* *********************************************************************************************************//**
 * Returns the real file name for the (a_index+1)^th Protare. The index **a_index** can only be 0 (normal protare) or 1 (TNSL protare).
 *
 * @param           a_index [in]    The index of the Protare whose real file name is returned.
 *
 * @return                          The real file name.
 ******************************************************************/

std::string const &ProtareTNSL::realFileName( std::size_t a_index ) const {

    if( a_index == 0 ) return( m_protare->realFileName( ) );
    if( a_index == 1 ) return( m_TNSL->realFileName( ) );
    throw Exception( "ProtareTNSL::realFileName: index can only be 0 or 1." );
}

/* *********************************************************************************************************//**
 * Returns the list of libraries for the (a_index+1)^th contained Protare. The index **a_index** can only be 0 (normal protare) or 1 (TNSL protare).
 * 
 * @param           a_index     [in]    The index of the Protare whose libraries are to be returned.
 *
 * @return                              The list of libraries.
 ******************************************************************/

std::vector<std::string> ProtareTNSL::libraries( std::size_t a_index ) const {

    if( a_index == 0 ) return( m_protare->libraries( ) );
    if( a_index == 1 ) return( m_TNSL->libraries( ) );
    throw Exception( "ProtareTNSL::libraries: index can only be 0 or 1." );
}

/* *********************************************************************************************************//**
 * Returns the evaluation for the (a_index+1)^th Protare. The index **a_index** can only be 0 (normal protare) or 1 (TNSL protare).
 *
 * @param           a_index [in]    The index of the Protare whose evaluation is returned.
 *
 * @return                          The evaluation.
 ******************************************************************/

std::string const &ProtareTNSL::evaluation( std::size_t a_index ) const {

    if( a_index == 0 ) return( m_protare->evaluation( ) );
    if( a_index == 1 ) return( m_TNSL->evaluation( ) );
    throw Exception( "ProtareTNSL::evaluation: index can only be 0 or 1." );
}

/* *********************************************************************************************************//**
 * Returns the projectile frame for the (a_index+1)^th Protare. The index **a_index** can only be 0 (normal protare) or 1 (TNSL protare).
 *
 * @param           a_index [in]    The index of the Protare whose projectile frame is returned.
 *
 * @return                          The projectile frame.
 ******************************************************************/

Frame ProtareTNSL::projectileFrame( std::size_t a_index ) const {

    if( a_index == 0 ) return( m_protare->projectileFrame( ) );
    if( a_index == 1 ) return( m_TNSL->projectileFrame( ) );
    throw Exception( "ProtareTNSL::projectileFrame: index can only be 0 or 1." );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the (a_index - 1)th **ProtareSingle**.
 *
 * @param a_index               [in]    Index of the **ProtareSingle** to return. Can only be 0 or 1.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

ProtareSingle *ProtareTNSL::protare( std::size_t a_index ) {

    if( a_index == 0 ) return( m_protare );
    if( a_index == 1 ) return( m_TNSL );
    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the pointer representing the (a_index - 1)th **ProtareSingle**.
 *
 * @param a_index               [in]    Index of the **ProtareSingle** to return. Can only be 0 or 1.
 *
 * @return                              Pointer to the requested protare or nullptr if invalid *a_index*..
 ***********************************************************************************************************/

ProtareSingle const *ProtareTNSL::protare( std::size_t a_index ) const {

    if( a_index == 0 ) return( m_protare );
    if( a_index == 1 ) return( m_TNSL );
    return( nullptr );
}

/* *********************************************************************************************************//**
 * Returns the number of LazyParsingHelperForms instantiated.
 *
 * @return              The number of LazyParsingHelperForms instantiated.
 ******************************************************************/

int ProtareTNSL::numberOfLazyParsingHelperForms( ) const {

    return( m_protare->numberOfLazyParsingHelperForms( ) + m_TNSL->numberOfLazyParsingHelperForms( ) );
}

/* *********************************************************************************************************//**
 * Returns the number of instantiated LazyParsingHelperForms replaced with the appropriate form.
 * 
 * @return              The number of LazyParsingHelperForms replaced.
 ******************************************************************/
 
int ProtareTNSL::numberOfLazyParsingHelperFormsReplaced( ) const {

    return( m_protare->numberOfLazyParsingHelperFormsReplaced( ) + m_TNSL->numberOfLazyParsingHelperFormsReplaced( ) );
}

/* *********************************************************************************************************//**
 * Returns the threshold factor for the projectile hitting the target.
 *
 * @return              The threshold factor.
 ******************************************************************/

double ProtareTNSL::thresholdFactor( ) const {

    return( m_protare->thresholdFactor( ) );
}

/* *********************************************************************************************************//**
 * Returns the Documentation_1_10::Suite from the non TNSL protare.
 *
 * @return              The Documentation_1_10::Suite.
 ******************************************************************/

Documentation_1_10::Suite &ProtareTNSL::documentations( ) {

    return( m_protare->documentations( ) );
}

/* *********************************************************************************************************//**
 * Returns the style with label *a_label* from the non TNSL protare.
 *
 * @param           a_label     [in]    The label of the requested style.
 * @return                              The style with label *a_label*.
 ******************************************************************/

Styles::Base &ProtareTNSL::style( std::string const &a_label ) {

    return( m_protare->style( a_label ) );
}

/* *********************************************************************************************************//**
 * Returns the Styles::Suite from the non TNSL protare.
 *
 * @return              The Styles::Suite.
 ******************************************************************/

Styles::Suite &ProtareTNSL::styles( ) {

    return( m_protare->styles( ) );
}

/* *********************************************************************************************************//**
 * Returns the Styles::Suite from the non TNSL protare.
 *
 * @return              The Styles::Suite.
 ******************************************************************/

Styles::Suite const &ProtareTNSL::styles( ) const {

    return( m_protare->styles( ) );
}

/* *********************************************************************************************************//**
 * Returns the intid for the requested particle or -1 if the particle is not in *m_protare* PoPs database.
 *
 * @param a_id                 [in]    The GNDS PoPs id for particle whose intd is requested.
 *
 * @return                             C++ int for the requested particle or -1 if particle is not in PoPs. 
 ******************************************************************/

int ProtareTNSL::intid( std::string const &a_id ) const {

    return( m_protare->intid( a_id ) );
}

/* *********************************************************************************************************//**
 * Calls productIDs for each Protare contained in *this*.
 *
 * @param   a_ids                   [in]        Contains the list of particle ids.
 * @param   a_particles             [in]        The list of particles to be transported.
 * @param   a_transportablesOnly    [in]        If *true* only transportable particle ids are added to *a_ids*.
 ******************************************************************/

void ProtareTNSL::productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const {

    m_protare->productIDs( a_ids, a_particles, a_transportablesOnly );
    m_TNSL->productIDs( a_ids, a_particles, a_transportablesOnly );
}

/* *********************************************************************************************************//**
 * Determines the maximum Legredre order present in the multi-group transfer matrix for a give product for a give label.
 * Loops over all contained Protares to determine the maximum Legredre order.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    The id of the requested product.
 *
 * @return                              The maximum Legredre order. If no transfer matrix data are present for the requested product, -1 is returned.
 ***********************************************************************************************************/

int ProtareTNSL::maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    int maximumLegendreOrder1 = m_protare->maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID );
    int maximumLegendreOrder2 = m_TNSL->maximumLegendreOrder( a_smr, a_settings, a_temperatureInfo, a_productID );

    if( maximumLegendreOrder1 > maximumLegendreOrder2 ) return( maximumLegendreOrder1 );
    return( maximumLegendreOrder2 );
}

/* *********************************************************************************************************//**
 * Returns a list of all process temperature data from the non TNSL protare. For each temeprature, the labels for its
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

Styles::TemperatureInfos ProtareTNSL::temperatures( ) const {

    return( m_protare->temperatures( ) );
}

/* *********************************************************************************************************//**
 * Returns the number of reactions from the non TNSL protare.
 *
 * @return              The total number of reactions.
 ******************************************************************/

std::size_t ProtareTNSL::numberOfReactions( ) const {

    return( m_TNSL->numberOfReactions( ) + m_protare->numberOfReactions( ) );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th reaction from the non TNSL protare.
 *
 * @param a_index           [in]    The index of the requested reaction.
 * @return                          The (*a_index*+1)th reaction.
 ***********************************************************************************************************/

Reaction *ProtareTNSL::reaction( std::size_t a_index ) {

    if( a_index < m_TNSL->numberOfReactions( ) ) return( m_TNSL->reaction( a_index ) );

    return( m_protare->reaction( a_index - m_TNSL->numberOfReactions( ) ) );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th reaction from the non TNSL protare.
 *      
 * @param a_index           [in]    The index of the requested reaction.
 * @return                          The (*a_index*+1)th reaction.
 ***********************************************************************************************************/
    
Reaction const *ProtareTNSL::reaction( std::size_t a_index ) const {

    if( a_index < m_TNSL->numberOfReactions( ) ) return( m_TNSL->reaction( a_index ) );

    return( m_protare->reaction( a_index - m_TNSL->numberOfReactions( ) ) );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th reaction from the non TNSL protare or **nullptr**. If the indexed reaction is 
 * deactivated or exlucded, a **nullptr** is returned.
 *
 * @param a_index               [in]    The index of the requested reaction.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_reactionsToExclude  [in]    A list of reaction indices that are to be ignored when calculating the cross section.
 *
 * @return                          The (*a_index*+1)th reaction or a **nullptr**.
 ***********************************************************************************************************/

Reaction const *ProtareTNSL::reaction( std::size_t a_index, Transporting::MG const &a_settings, 
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    if( a_index < m_TNSL->numberOfReactions( ) ) return( m_TNSL->reaction( a_index, a_settings, a_reactionsToExclude ) );

    return( m_protare->reaction( a_index - m_TNSL->numberOfReactions( ), a_settings, a_reactionsToExclude ) );
}

/* *********************************************************************************************************//**
 * Returns the number of orphanProduct's from the non TNSL protare.
 *  
 * @return              The total number of orphanProducts.
 ******************************************************************/
    
std::size_t ProtareTNSL::numberOfOrphanProducts( ) const {

    return( m_protare->numberOfOrphanProducts( ) );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th orphanProduct from the non TNSL protare.
 *
 * @param a_index           [in]    The index of the requested orphanProduct.
 * @return                          The (*a_index*+1)th orphanProduct.
 ***********************************************************************************************************/

Reaction *ProtareTNSL::orphanProduct( std::size_t a_index ) {

    return( m_protare->orphanProduct( a_index ) );
}

/* *********************************************************************************************************//**
 * Returns the (*a_index*+1)th orphanProduct from the non TNSL protare.
 *
 * @param a_index           [in]    The index of the requested orphanProduct.
 * @return                          The (*a_index*+1)th orphanProduct.
 ***********************************************************************************************************/

Reaction const *ProtareTNSL::orphanProduct( std::size_t a_index ) const {

    return( m_protare->orphanProduct( a_index ) );
}

/* *********************************************************************************************************//**
 * Re-indexs the reactions in the reactions, orphanProducts and fissionComponents suites.
 *
 ***********************************************************************************************************/
    
void ProtareTNSL::updateReactionIndices( LUPI_maybeUnused std::size_t a_offset ) const {

    m_TNSL->updateReactionIndices( 0 );
    m_protare->updateReactionIndices( m_TNSL->numberOfReactions( ) );
}

/* *********************************************************************************************************//**
 * Returns true if at least one reaction contains a fission channel.
 *
 * @return      true if at least one reaction contains a fission channel and false otherwise.
 ***********************************************************************************************************/

bool ProtareTNSL::hasFission( ) const {

    return( m_protare->hasFission( ) );
}

/* *********************************************************************************************************//**
 * Returns **false* if protare has delayed fission neutrons for an active reaction and they are not complete; otherwise, returns **true**.
 *  
 * @return      bool        
 ***********************************************************************************************************/
 
bool ProtareTNSL::isDelayedFissionNeutronComplete( ) const {
 
    return( m_protare->isDelayedFissionNeutronComplete( ) );
}

/* *********************************************************************************************************//**
 * Returns the multi-group boundaries for the requested label and product.
 *
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 * @param a_productID           [in]    ID for the requested product.
 *
 * @return                              List of multi-group boundaries.
 ***********************************************************************************************************/

std::vector<double> ProtareTNSL::groupBoundaries( Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const {

    return( m_protare->groupBoundaries( a_settings, a_temperatureInfo, a_productID ) );
}

/* *********************************************************************************************************//**
 * Returns the inverse speeds for the requested label from the non TNSL protare. The label must be for a heated multi-group style.
 *
 * @param a_smr                 [Out]   If errors are not to be thrown, then the error is reported via this instance.
 * @param a_settings            [in]    Specifies the requested label.
 * @param a_temperatureInfo     [in]    Specifies the temperature and labels use to lookup the requested data.
 *
 * @return                              List of inverse speeds.
 ***********************************************************************************************************/

Vector ProtareTNSL::multiGroupInverseSpeed( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo ) const {

    return( m_protare->multiGroupInverseSpeed( a_smr, a_settings, a_temperatureInfo ) );
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

Vector ProtareTNSL::multiGroupCrossSection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude,
                std::string const &a_label ) const {

    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Vector vector = m_protare->multiGroupCrossSection( a_smr, a_settings, a_temperatureInfo, excludeReactionsSet, a_label );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( vector );

    Vector vectorElastic = m_elasticReaction->multiGroupCrossSection( a_smr, a_settings, a_temperatureInfo, a_label );

    combineVectors( a_settings, a_temperatureInfo, vector, vectorElastic, 
            m_TNSL->multiGroupCrossSection( a_smr, a_settings, a_temperatureInfo, excludeReactionsSet, a_label ) );
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

Vector ProtareTNSL::multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, bool a_final, bool a_effectivePhotoAtomic,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Vector vector = m_protare->multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final, a_effectivePhotoAtomic );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( vector );

    Vector vectorElastic = m_elasticReaction->multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final );

    combineVectors( a_settings, a_temperatureInfo, vector, vectorElastic, m_TNSL->multiGroupQ( a_smr, a_settings, a_temperatureInfo, a_final, a_effectivePhotoAtomic ) );
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

Vector ProtareTNSL::multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Vector vector = m_protare->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( vector );

    Vector vectorElastic = m_elasticReaction->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID );

    combineVectors( a_settings, a_temperatureInfo, vector, vectorElastic, m_TNSL->multiGroupMultiplicity( a_smr, a_settings, a_temperatureInfo, a_productID ) );
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

Vector ProtareTNSL::multiGroupFissionNeutronMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    return( m_protare->multiGroupFissionNeutronMultiplicity( a_smr, a_settings, a_temperatureInfo, a_reactionsToExclude ) );
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

Vector ProtareTNSL::multiGroupFissionGammaMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    return( m_protare->multiGroupFissionGammaMultiplicity( a_smr, a_settings, a_temperatureInfo, a_reactionsToExclude ) );
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

Matrix ProtareTNSL::multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, 
                std::string const &a_productID, std::size_t a_order, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Matrix matrix = m_protare->multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order, excludeReactionsSet );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( matrix );

    Matrix matrixElastic = m_elasticReaction->multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order );
    Matrix matrixTNSL = m_TNSL->multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_productID, a_order );

    combineMatrices( a_settings, a_temperatureInfo, matrix, matrixElastic, matrixTNSL );
    return( matrix );
}

/* *********************************************************************************************************//**
 * Like ProtareTNSL::multiGroupProductMatrix, but only returns the fission neutron, transfer matrix.
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

Matrix ProtareTNSL::multiGroupFissionMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::size_t a_order,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    return( m_protare->multiGroupFissionMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, a_order, a_reactionsToExclude ) );
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

Vector ProtareTNSL::multiGroupTransportCorrection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::size_t a_order, 
                TransportCorrectionType a_transportCorrectionType, double a_temperature, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    if( a_transportCorrectionType == TransportCorrectionType::None ) return( Vector( 0 ) );

    Matrix matrix( multiGroupProductMatrix( a_smr, a_settings, a_temperatureInfo, a_particles, projectile( ).ID( ), a_order + 1, 
            a_reactionsToExclude ) );
    Matrix matrixCollapsed = collapse( matrix, a_settings, a_particles, a_temperature, projectile( ).ID( ) );
    std::size_t size = matrixCollapsed.size( );
    std::vector<double> transportCorrection1( size, 0 );

    if( a_transportCorrectionType == TransportCorrectionType::None ) {
        }
    else if( a_transportCorrectionType == TransportCorrectionType::Pendlebury ) {
        for( std::size_t index = 0; index < size; ++index ) transportCorrection1[index] = matrixCollapsed[index][index]; }
    else {
        throw Exception( "Unsupported transport correction: only None and Pendlebury (i.e., Pendlebury/Underhill) are currently supported." );
    }
    return( Vector( transportCorrection1 ) );
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

Vector ProtareTNSL::multiGroupAvailableEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Vector vector = m_protare->multiGroupAvailableEnergy( a_smr, a_settings, a_temperatureInfo, excludeReactionsSet );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( vector );

    Vector vectorElastic = m_elasticReaction->multiGroupAvailableEnergy( a_smr, a_settings, a_temperatureInfo );

    combineVectors( a_settings, a_temperatureInfo, vector, vectorElastic, m_TNSL->multiGroupAvailableEnergy( a_smr, a_settings, a_temperatureInfo ) );
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

Vector ProtareTNSL::multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Vector vector = m_protare->multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID, excludeReactionsSet );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( vector );

    Vector vectorElastic = m_elasticReaction->multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID );

    combineVectors( a_settings, a_temperatureInfo, vector, vectorElastic, m_TNSL->multiGroupAverageEnergy( a_smr, a_settings, a_temperatureInfo, a_productID ) );
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

Vector ProtareTNSL::multiGroupDepositionEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Vector vector = m_protare->multiGroupDepositionEnergy( a_smr, a_settings, a_temperatureInfo, a_particles, excludeReactionsSet );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( vector );

    Vector vectorElastic = m_elasticReaction->multiGroupDepositionEnergy( a_smr, a_settings, a_temperatureInfo, a_particles );

    combineVectors( a_settings, a_temperatureInfo, vector, vectorElastic, m_TNSL->multiGroupDepositionEnergy( a_smr, a_settings, a_temperatureInfo, a_particles ) );
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

Vector ProtareTNSL::multiGroupAvailableMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude ) const {

    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Vector vector = m_protare->multiGroupAvailableMomentum( a_smr, a_settings, a_temperatureInfo, excludeReactionsSet );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( vector );

    Vector vectorElastic = m_elasticReaction->multiGroupAvailableMomentum( a_smr, a_settings, a_temperatureInfo );

    combineVectors( a_settings, a_temperatureInfo, vector, vectorElastic, m_TNSL->multiGroupAvailableMomentum( a_smr, a_settings, a_temperatureInfo ) );
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

Vector ProtareTNSL::multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Vector vector = m_protare->multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID, excludeReactionsSet );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( vector );

    Vector vectorElastic = m_elasticReaction->multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID );

    combineVectors( a_settings, a_temperatureInfo, vector, vectorElastic, m_TNSL->multiGroupAverageMomentum( a_smr, a_settings, a_temperatureInfo, a_productID ) );
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

Vector ProtareTNSL::multiGroupDepositionMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Vector vector = m_protare->multiGroupDepositionMomentum( a_smr, a_settings, a_temperatureInfo, a_particles, excludeReactionsSet );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( vector );

    Vector vectorElastic = m_elasticReaction->multiGroupDepositionMomentum( a_smr, a_settings, a_temperatureInfo, a_particles );

    combineVectors( a_settings, a_temperatureInfo, vector, vectorElastic, m_TNSL->multiGroupDepositionMomentum( a_smr, a_settings, a_temperatureInfo, a_particles ) );
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

Vector ProtareTNSL::multiGroupGain( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                ExcludeReactionsSet const &a_reactionsToExclude ) const {

    std::string const projectile_id = m_protare->projectile( ).ID( );    
    ExcludeReactionsSet excludeReactionsSet( a_reactionsToExclude );
    Vector vector = m_protare->multiGroupGain( a_smr, a_settings, a_temperatureInfo, a_productID, excludeReactionsSet );
    excludeReactionsSetAdjust( excludeReactionsSet, *m_protare );

    if( !m_elasticReaction->active( ) ) return( vector );

    Vector vectorElastic = m_elasticReaction->multiGroupGain( a_smr, a_settings, a_temperatureInfo, a_productID, projectile_id );

    combineVectors( a_settings, a_temperatureInfo, vector, vectorElastic, m_TNSL->multiGroupGain( a_smr, a_settings, a_temperatureInfo, a_productID ) );

    return( vector );
}

/* *********************************************************************************************************//**
 * If the protare is a ProtareTNSL then summing over all reactions will include the standard protare's elastic cross section
 * in the domain of the TNSL data. The standard elastic cross section should not be added in this domain.
 * If needed, this function corrects the cross section for this over counting of the elastic cross section.
 *
 * @param       a_label                     [in]    The label of the elastic cross section data to use if over counting needs to be corrected.
 * @param       a_crossSectionSum           [in]    The cross section to correct.
 ***********************************************************************************************************/

void ProtareTNSL::TNSL_crossSectionSumCorrection( std::string const &a_label, Functions::XYs1d &a_crossSectionSum ) {

    double projectileEnergyMax = m_TNSL->projectileEnergyMax( );
    Functions::XYs1d *xys1d = m_elasticReaction->crossSection( ).get<Functions::XYs1d>( a_label );

    ptwXYPoints *ptwXY = const_cast<ptwXYPoints *>( xys1d->ptwXY( ) );
    ptwXYPoints *sliced = ptwXY_domainMaxSlice( NULL, ptwXY, projectileEnergyMax, ptwXY_length( NULL, ptwXY ), 1 );
    Functions::XYs1d slicedXYs1d( a_crossSectionSum.axes( ), sliced );

    a_crossSectionSum -= slicedXYs1d;
}

/* *********************************************************************************************************//**
 * If the protare is a ProtareTNSL then summing over all reactions will include the standard protare's elastic cross section
 * in the domain of the TNSL data. The standard elastic cross section should not be added in this domain.
 * If needed, this function corrects the cross section for this over counting of the elastic cross section.
 *
 * @param       a_label                     [in]    The label of the elastic cross section data to use if over counting needs to be corrected.
 * @param       a_crossSectionSum           [in]    The cross section to correct.
 ***********************************************************************************************************/

void ProtareTNSL::TNSL_crossSectionSumCorrection( std::string const &a_label, Functions::Ys1d &a_crossSectionSum ) {

    double projectileEnergyMax = m_TNSL->projectileEnergyMax( );
    Styles::GriddedCrossSection const *griddedCrossSection = m_protare->styles( ).get<Styles::GriddedCrossSection>( a_label );
    nf_Buffer<double> const energies = griddedCrossSection->grid( ).values( );
    Functions::Ys1d const *ys1d = m_elasticReaction->crossSection( ).get<Functions::Ys1d const>( a_label );

    std::size_t start = ys1d->start( );
    std::vector<double> const &Ys = ys1d->Ys( );
    std::vector<double> &crossSectionYs = a_crossSectionSum.Ys( );
    for( std::size_t index = 0; index < Ys.size( ); ++index ) {
        if( energies[index] > projectileEnergyMax ) break;
        crossSectionYs[index+start] -= Ys[index];
    }
}

/* *********************************************************************************************************//**
 * This method always returns 1 since the projectile is always a neutron.
 *
 * @return      Always returns 1.
 ***********************************************************************************************************/

stringAndDoublePairs ProtareTNSL::muCutoffForCoulombPlusNuclearElastic( ) const {

    stringAndDoublePairs stringAndDoublePairs1;

    return( stringAndDoublePairs1 );
}

/* *********************************************************************************************************//**
 * Calls the **incompleteParticles** method for each **ProtareSingle** in *this*.
 *
 * @param a_settings                [in]    Specifies the requested label.
 * @param a_incompleteParticles     [out]   The list of particles whose **completeParticle** method returns *false*.
 ***********************************************************************************************************/
 
void ProtareTNSL::incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const {
    
    m_protare->incompleteParticles( a_settings, a_incompleteParticles );
    m_TNSL->incompleteParticles( a_settings, a_incompleteParticles );
}

}
