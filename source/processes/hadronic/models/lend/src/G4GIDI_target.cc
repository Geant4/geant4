//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include <G4GIDI.hh>

/*! \class G4GIDI_target
 * A class to store map files for a particular projectile.
 */

/* *********************************************************************************************************//**
 * @param a_MCProtare               [in]    The MCGIDI protare.
 ***********************************************************************************************************/

G4GIDI_target::G4GIDI_target( PoPI::Database const &a_pops, MCGIDI::DomainHash const &a_domainHash, GIDI::Protare const &a_GIDI_protare, 
                MCGIDI::Protare *a_MCGIDI_protare ) :
        m_MCGIDI_protare( a_MCGIDI_protare ),
        m_target( a_GIDI_protare.target( ).ID( ) ),
        m_fileName( a_GIDI_protare.fileName( ) ),
        m_targetZ( 0 ),
        m_targetA( 0 ),
        m_targetM( 0 ),
        m_targetMass( 0.0 ),
        m_domainHash( a_domainHash ),
        m_elasticAngular( nullptr ) {

    PoPI::Base const *targetAsBase = &a_pops.get<PoPI::Base const>( m_target );
    PoPI::Base const *targetAsBase2 = targetAsBase;
    if( targetAsBase->isAlias( ) ) {
        targetAsBase2 = &a_pops.get<PoPI::Base const>( a_pops.final( m_target ) );
    }
    m_targetZ = PoPI::particleZ( *targetAsBase2, true );
    m_targetA = PoPI::particleA( *targetAsBase2, true );
    m_targetM = PoPI::particleMetaStableIndex( *targetAsBase );
    if( targetAsBase2->isParticle( ) ) {
        PoPI::Particle const *targetAsParticle = static_cast<PoPI::Particle const *>( targetAsBase2 );
        m_targetMass = targetAsParticle->massValue( "amu" );
    }

    for( std::size_t reactionIndex = 0; reactionIndex < a_MCGIDI_protare->numberOfReactions( ); ++reactionIndex ) {
        MCGIDI::Reaction const *reaction = a_MCGIDI_protare->reaction( reactionIndex );

        if( reaction->ENDF_MT( ) == 2 ) {
            m_elasticIndices.push_back( static_cast<int>( reactionIndex ) ); }
        else if( reaction->ENDF_MT( ) == 18 ) {
            m_fissionIndices.push_back( static_cast<int>( reactionIndex ) ); }
        else if( reaction->ENDF_MT( ) == 102 ) {
            m_captureIndices.push_back( static_cast<int>( reactionIndex ) ); }
        else {
            m_othersIndices.push_back( static_cast<int>( reactionIndex ) );
        }
    }

    if( m_elasticIndices.size( ) > 0 ) {
        MCGIDI::Reaction const *elastic = a_MCGIDI_protare->reaction( m_elasticIndices[0] );
        MCGIDI::Product const *firstProduct = elastic->product( 0 );
        MCGIDI::Distributions::AngularTwoBody const *angularTwoBody = static_cast<MCGIDI::Distributions::AngularTwoBody const *>( firstProduct->distribution( ) );
        m_elasticAngular = angularTwoBody->angular( );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

G4GIDI_target::~G4GIDI_target( ) {

    delete m_MCGIDI_protare;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

int G4GIDI_target::getNumberOfChannels( ) const {

    return( static_cast<int>( m_MCGIDI_protare->numberOfReactions( ) ) );
}
/* *********************************************************************************************************//**
 ***********************************************************************************************************/

int G4GIDI_target::getNumberOfProductionChannels( ) const {

    return( static_cast<int>( m_MCGIDI_protare->numberOfOrphanProducts( ) ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

channelID G4GIDI_target::getChannelsID( int channelIndex ) const {

    return( channelID( m_MCGIDI_protare->reaction( channelIndex )->label( ).c_str( ) ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

std::vector<channelID> *G4GIDI_target::getChannelIDs( ) const {

    std::vector<channelID> *channelIDs = new std::vector<channelID>( 0 );

    for( std::size_t reactionIndex = 0; reactionIndex < m_MCGIDI_protare->numberOfReactions( ); ++reactionIndex ) {
        MCGIDI::Reaction const *reaction = m_MCGIDI_protare->reaction( reactionIndex );
        channelIDs->push_back( reaction->label( ).c_str( ) );
    }

    return( channelIDs );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

std::vector<channelID> *G4GIDI_target::getProductionChannelIDs( ) const {

    std::vector<channelID> *channelIDs = new std::vector<channelID>( 0 );

    for( std::size_t reactionIndex = 0; reactionIndex < m_MCGIDI_protare->numberOfOrphanProducts( ); ++reactionIndex ) {
        MCGIDI::Reaction const *reaction = m_MCGIDI_protare->orphanProduct( reactionIndex );
        channelIDs->push_back( reaction->label( ).c_str( ) );
    }

    return( channelIDs );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

double G4GIDI_target::getTotalCrossSectionAtE( double a_energy, double a_temperature ) const {

    int hashIndex = m_domainHash.index( a_energy );

    return( m_MCGIDI_protare->crossSection( m_URR_protareInfos, hashIndex, a_temperature, a_energy ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

double G4GIDI_target::getElasticCrossSectionAtE( double a_energy, double a_temperature ) const {

    return( sumChannelCrossSectionAtE( m_elasticIndices, a_energy, a_temperature ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

double G4GIDI_target::getCaptureCrossSectionAtE( double a_energy, double a_temperature ) const {

    return( sumChannelCrossSectionAtE( m_captureIndices, a_energy, a_temperature ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

double G4GIDI_target::getFissionCrossSectionAtE( double a_energy, double a_temperature ) const {

    return( sumChannelCrossSectionAtE( m_fissionIndices, a_energy, a_temperature ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

double G4GIDI_target::getOthersCrossSectionAtE( double a_energy, double a_temperature ) const {

    return( sumChannelCrossSectionAtE( m_othersIndices, a_energy, a_temperature ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

double G4GIDI_target::sumChannelCrossSectionAtE( std::vector<int> const &a_indices, double a_energy, double a_temperature ) const {

    int hashIndex = m_domainHash.index( a_energy );
    double crossSection = 0.0;

    for( auto indexIter = a_indices.begin( ); indexIter != a_indices.end( ); ++indexIter ) {
        crossSection += m_MCGIDI_protare->reactionCrossSection( *indexIter, m_URR_protareInfos, hashIndex, a_temperature, a_energy );
    }

    return( crossSection );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

double G4GIDI_target::sumChannelCrossSectionAtE( int a_nIndices, int const *a_indices, double a_energy, double a_temperature ) const {

    std::vector<int> indices( a_nIndices );
    for( int index = 0; index < a_nIndices; ++index ) indices[index] = a_indices[index];

    return( sumChannelCrossSectionAtE( indices, a_energy, a_temperature ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

int G4GIDI_target::sampleChannelCrossSectionAtE( std::vector<int> const &a_indices, double a_energy, double a_temperature, 
                double (*a_rng)( void * ), void *a_rngState ) const {

    int index = 0;
    int nIndices = static_cast<int>( a_indices.size( ) );
    double crossSectionSample = sumChannelCrossSectionAtE( a_indices, a_energy, a_temperature );
    double crossSectionSum = 0.0;

    crossSectionSample *= a_rng( a_rngState );
    for( ; index < nIndices; ++index ) {
        crossSectionSum += sumChannelCrossSectionAtE( 1, &a_indices[index], a_energy, a_temperature );
        if( crossSectionSum >= crossSectionSample ) break;
    }
    if( index == nIndices ) --index;          // This should really be an error.
    
    return( a_indices[index] );
}
/* *********************************************************************************************************//**
 ***********************************************************************************************************/

int G4GIDI_target::sampleChannelCrossSectionAtE( int a_nIndices, int const *a_indices, double a_energy, double a_temperature, 
                double (*a_rng)( void * ), void *a_rngState ) const {

    std::vector<int> indices( a_nIndices );
    for( int index = 0; index < a_nIndices; ++index ) indices[index] = a_indices[index];

    return( sampleChannelCrossSectionAtE( indices, a_energy, a_temperature, a_rng, a_rngState ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

double G4GIDI_target::getElasticFinalState( double a_energy, LUPI_maybeUnused double a_temperature, double (*a_rng)( void * ), void *a_rngState ) const {

    return( m_elasticAngular->sample( a_energy, a_rng( a_rngState ), [&]() -> double { return a_rng( a_rngState ); } ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

std::vector<G4GIDI_Product> *G4GIDI_target::getCaptureFinalState( double a_energy, double a_temperature, double (*a_rng)( void * ), void *a_rngState ) const {

    return( getFinalState( m_captureIndices, a_energy, a_temperature, a_rng, a_rngState ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

std::vector<G4GIDI_Product> *G4GIDI_target::getFissionFinalState( double a_energy, double a_temperature, double (*a_rng)( void * ), void *a_rngState ) const {

    return( getFinalState( m_fissionIndices, a_energy, a_temperature, a_rng, a_rngState ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

std::vector<G4GIDI_Product> *G4GIDI_target::getOthersFinalState( double a_energy, double a_temperature, double (*a_rng)( void * ), void *a_rngState ) const {

    return( getFinalState( m_othersIndices, a_energy, a_temperature, a_rng, a_rngState ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

std::vector<G4GIDI_Product> *G4GIDI_target::getFinalState( std::vector<int> const &a_indices, double a_energy, double a_temperature,
                        double (*a_rng)( void * ), void *a_rngState ) const {

    int reactionIndex;

    if( a_indices.size( ) == 0 ) return( NULL );
    if( a_indices.size( ) == 1 ) {
        reactionIndex = a_indices[0]; }
    else {
        reactionIndex = sampleChannelCrossSectionAtE( a_indices, a_energy, a_temperature, a_rng, a_rngState );
    }

    MCGIDI::Sampling::StdVectorProductHandler productHandler;
    MCGIDI::Sampling::Input input( false, MCGIDI::Sampling::Upscatter::Model::none );

    MCGIDI::Reaction const *reaction = m_MCGIDI_protare->reaction( reactionIndex );
    reaction->sampleProducts( m_MCGIDI_protare, a_energy, input, [&]() -> double { return a_rng( a_rngState ); },
        [&] (MCGIDI::Sampling::Product &a_product) -> void { productHandler.push_back( a_product ); }, productHandler );

    std::vector<G4GIDI_Product> *products = new std::vector<G4GIDI_Product>( productHandler.size( ) );

    for( std::size_t index = 0; index < productHandler.size( ); ++index ) {
        MCGIDI::Sampling::Product &productIn = productHandler[index];
        G4GIDI_Product &productOut = (*products)[index];

        if( productIn.m_productIntid  == 1020000000 ) {             // Neutron.
            productOut.A = 1;
            productOut.Z = 0;
            productOut.m = 0; }
        else if( productIn.m_productIntid  == 1000000000 ) {        // Photon.
            productOut.A = 0;
            productOut.Z = 0;
            productOut.m = 0; }
        else {
            PoPI::ParseIntidInfo parseIntidInfo( productIn.m_productIntid );
            productOut.A = parseIntidInfo.AAA( );
            productOut.Z = parseIntidInfo.ZZZ( );
            productOut.m = parseIntidInfo.metaStableIndex( );
        }

        productOut.kineticEnergy = productIn.m_kineticEnergy;
        productOut.px = productIn.m_px_vx;
        productOut.py = productIn.m_py_vy;
        productOut.pz = productIn.m_pz_vz;
    }

    return( products );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

std::vector<G4GIDI_Product> *G4GIDI_target::getFinalState( int a_nIndices, int const *a_indices, double a_energy, double a_temperature, 
                        double (*a_rng)( void * ), void *a_rngState ) const {

    std::vector<int> indices( a_nIndices );
    for( int index = 0; index < a_nIndices; ++index ) indices[index] = a_indices[index];

    return( getFinalState( indices, a_energy, a_temperature, a_rng, a_rngState ) );
}
