/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <iostream>
#include <stdlib.h>

#include "GIDI_settings.hh"

/*
=========================================================
*/
/**
    This is the top settings class used when a GND file is read.
*/
GIDI_settings::GIDI_settings( ) {

}
/*
=========================================================
*/
GIDI_settings::~GIDI_settings( ) {

}
/*
=========================================================
*/
int GIDI_settings::addParticle( GIDI_settings_particle const &particle ) {

    int PoPId = particle.getPoPId( );

    if( mParticles.find( PoPId ) != mParticles.end( ) ) return( 1 );
    mParticles.insert( std::pair<int, GIDI_settings_particle>( PoPId, GIDI_settings_particle( particle ) ) );
    return( 0 );
}
/*
=========================================================
*/
GIDI_settings_particle const *GIDI_settings::getParticle( int PoPId ) const {

    std::map<int, GIDI_settings_particle>::const_iterator particle = mParticles.find( PoPId );

    if( particle == mParticles.end( ) ) return( NULL );
    return( &(particle->second) );
}
/*
=========================================================
*/
int GIDI_settings::eraseParticle( int PoPId ) {

    std::map<int, GIDI_settings_particle>::iterator particle = mParticles.find( PoPId );

    if( particle == mParticles.end( ) ) return( 1 );
    mParticles.erase( PoPId );
    return( 0 );
}
