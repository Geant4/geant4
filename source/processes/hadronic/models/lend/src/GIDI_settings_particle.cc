/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <iostream>
#include <cmath>

#include "GIDI_settings.hh"

using namespace GIDI;

/*
=========================================================
*/
GIDI_settings_particle::GIDI_settings_particle( int PoPId, bool transporting, int energyMode ) : mGroup( ) {

    initialize( PoPId, transporting, energyMode );
}
/*
=========================================================
*/
GIDI_settings_particle::GIDI_settings_particle( GIDI_settings_particle const &particle ) {

    initialize( particle.mPoPId, particle.mTransporting, particle.mEnergyMode );
    setGroup( particle.mGroup );
    for( std::vector<GIDI_settings_processedFlux>::const_iterator iter = particle.mProcessedFluxes.begin( ); iter != particle.mProcessedFluxes.end( ); ++iter ) {
        mProcessedFluxes.push_back( *iter );
    }
}
/*
=========================================================
*/
int GIDI_settings_particle::initialize( int PoPId, bool transporting, int energyMode ) {

    mPoPId = PoPId;
    mTransporting = transporting;
    int energyMode_ =   ( energyMode & GIDI_settings_projectileEnergyMode_continuousEnergy )
                      + ( energyMode & GIDI_settings_projectileEnergyMode_grouped )
//                    + ( energyMode & GIDI_settings_projectileEnergyMode_fixedGrid ) // Currently not supported.
                        ;

    if( energyMode_ != energyMode ) throw 1;
    mEnergyMode = energyMode;

    mGroupX = NULL;
    setGroup( mGroup );
    return( 0 );
}
/*
=========================================================
*/
void GIDI_settings_particle::setGroup( GIDI_settings_group const &group ) {

    nfu_status status_nf;

    mGroup = group;

    if( mGroupX != NULL ) ptwX_free( mGroupX );
    mGroupX = NULL;
    if( mGroup.size( ) > 0 ) {
        if( ( mGroupX = ptwX_create( (int) mGroup.size( ), (int) mGroup.size( ), mGroup.pointer( ), &status_nf ) ) == NULL ) throw 1;
    }
}
/*
=========================================================
*/
GIDI_settings_particle::~GIDI_settings_particle( ) {

    if( mGroupX != NULL ) ptwX_free( mGroupX );
}
/*
=========================================================
*/
int GIDI_settings_particle::addFlux( statusMessageReporting* /* smr */, GIDI_settings_flux const &flux ) {

    double temperature = flux.getTemperature( );
    std::vector<GIDI_settings_processedFlux>::iterator iter;

    for( iter = mProcessedFluxes.begin( ); iter != mProcessedFluxes.end( ); ++iter ) {
        if( temperature <= iter->getTemperature( ) ) break;
    }
// BRB need to check if temperature is the same.
    mProcessedFluxes.insert( iter, GIDI_settings_processedFlux( flux, mGroupX ) );
    return( 0 );
}
/*
=========================================================
*/
GIDI_settings_processedFlux const *GIDI_settings_particle::nearestFluxToTemperature( double temperature ) const {

    double priorTemperature, lastTemperature;
    std::vector<GIDI_settings_processedFlux>::const_iterator iter;

    if( mProcessedFluxes.size( ) == 0 ) return( NULL );

    priorTemperature = mProcessedFluxes[0].getTemperature( );
    //TK adds next line
    lastTemperature  = mProcessedFluxes[0].getTemperature( );
    for( iter = mProcessedFluxes.begin( ); iter != mProcessedFluxes.end( ); ++iter ) {
        lastTemperature = iter->getTemperature( );
        if( lastTemperature > temperature ) break;
        //TK add next line
        priorTemperature = iter->getTemperature( );
    }
    if( iter == mProcessedFluxes.end( ) ) {
        --iter; }
    else {
        //if( fabs( lastTemperature - temperature ) < fabs( temperature - priorTemperature ) ) --iter;
        //TK modified above line
        if( std::fabs( lastTemperature - temperature ) > std::fabs( temperature - priorTemperature ) ) --iter;
    }
    return( &(*iter) );
}
/*
=========================================================
*/
ptwXPoints *GIDI_settings_particle::groupFunction( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double temperature, int order ) const {

    if( mGroupX == NULL ) return( NULL );
    GIDI_settings_processedFlux const *processedFlux = nearestFluxToTemperature( temperature );
    if( processedFlux == NULL ) return( NULL );
    return( processedFlux->groupFunction( smr, mGroupX, ptwXY1, order ) );
}

/*  ---- GIDI_settings_processedFlux ----  */
/*
=========================================================
*/
GIDI_settings_processedFlux::GIDI_settings_processedFlux( GIDI_settings_flux const &flux, ptwXPoints *groupX ) : mFlux( flux ) {

    nfu_status status_nf;
    ptwXYPoints *fluxXY = NULL;
    ptwXPoints *groupedFluxX;
    GIDI_settings_flux_order const *fluxOrder;
    double const *energies, *fluxes;

    for( int order = 0; order < (int) flux.size( ); ++order ) {
        fluxOrder = flux[order];
        energies = fluxOrder->getEnergies( );
        fluxes = fluxOrder->getFluxes( );
        if( ( fluxXY = ptwXY_createFrom_Xs_Ys( ptwXY_interpolationLinLin, NULL, 12, 1e-3, fluxOrder->size( ), 10, 
            fluxOrder->size( ), energies, fluxes, &status_nf, 0 ) ) == NULL ) goto err;
        mFluxXY.push_back( fluxXY );
        if( ( groupedFluxX = ptwXY_groupOneFunction( fluxXY, groupX, ptwXY_group_normType_none, NULL, &status_nf ) ) == NULL ) goto err;
        mGroupedFlux.push_back( groupedFluxX );
    }
    return;

err:
    throw 1;
}
/*
=========================================================
*/
GIDI_settings_processedFlux::GIDI_settings_processedFlux( GIDI_settings_processedFlux const &flux ) : mFlux( flux.mFlux ) {

    nfu_status status_nf;
    ptwXYPoints *fluxXY;
    ptwXPoints *fluxX;

    for( int order = 0; order < (int) mFlux.size( ); ++order ) {
        if( ( fluxXY = ptwXY_clone( flux.mFluxXY[order], &status_nf ) ) == NULL ) goto err;
        mFluxXY.push_back( fluxXY );
        if( ( fluxX = ptwX_clone( flux.mGroupedFlux[order], &status_nf ) ) == NULL ) goto err;
        mGroupedFlux.push_back( fluxX );
    }
    return;

err:
    for( std::vector<ptwXYPoints *>::iterator iter = mFluxXY.begin( ); iter != mFluxXY.end( ); ++iter ) ptwXY_free( *iter );
    for( std::vector<ptwXPoints *>::iterator iter = mGroupedFlux.begin( ); iter != mGroupedFlux.end( ); ++iter ) ptwX_free( *iter );
    throw 1;
}
/*
=========================================================
*/
GIDI_settings_processedFlux::~GIDI_settings_processedFlux( ) {

    for( std::vector<ptwXYPoints *>::iterator iter = mFluxXY.begin( ); iter != mFluxXY.end( ); ++iter ) ptwXY_free( *iter );
    for( std::vector<ptwXPoints *>::iterator iter = mGroupedFlux.begin( ); iter != mGroupedFlux.end( ); ++iter ) ptwX_free( *iter );
}
/*
=========================================================
*/
ptwXPoints *GIDI_settings_processedFlux::groupFunction( statusMessageReporting * /*smr*/, ptwXPoints *groupX, ptwXYPoints *ptwXY1, int order ) const {

    nfu_status status_nf;
    ptwXYPoints *fluxXY;
    ptwXPoints *groupedX;

    if( groupX == NULL ) return( NULL );
    if( order < 0 ) order = 0;
    if( order >= (int) mFluxXY.size( ) ) order = (int) mFluxXY.size( ) - 1;

    fluxXY = ptwXY_xSlice( mFluxXY[order], ptwXY_getXMin( ptwXY1 ), ptwXY_getXMax( ptwXY1 ), 10, 1, &status_nf );
//    if( fluxXY == NULL ) smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_xSlice error %d: %s", status_nf, nfu_statusMessage( status_nf ) );
    
    groupedX = ptwXY_groupTwoFunctions( ptwXY1, fluxXY, groupX, ptwXY_group_normType_norm, mGroupedFlux[order], &status_nf );
    ptwXY_free( fluxXY );
//    if( groupedX == NULL ) smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_groupTwoFunctions error %d: %s", status_nf, nfu_statusMessage( status_nf ) );
    return( groupedX );
}
