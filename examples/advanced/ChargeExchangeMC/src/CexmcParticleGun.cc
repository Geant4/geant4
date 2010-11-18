/*
 * ============================================================================
 *
 *       Filename:  CexmcParticleGun.cc
 *
 *    Description:  particle gun
 *
 *        Version:  1.0
 *        Created:  15.12.2009 13:29:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include "CexmcParticleGun.hh"
#include "CexmcParticleGunMessenger.hh"


CexmcParticleGun::CexmcParticleGun( CexmcPhysicsManager *  physicsManager,
                                    G4int  nmbOfParticles ) :
    G4ParticleGun( nmbOfParticles ), physicsManager( physicsManager ),
    origPos( 0, 0, 0 ), origDir( 0, 0, 0 ), origMomentumAmp( 0 ),
    messenger( NULL )
{
    messenger = new CexmcParticleGunMessenger( this );
}


CexmcParticleGun::~CexmcParticleGun()
{
    delete messenger;
}

