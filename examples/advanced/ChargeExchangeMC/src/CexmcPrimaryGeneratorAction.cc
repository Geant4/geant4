/*
 * =============================================================================
 *
 *       Filename:  CexmcPrimaryGeneratorAction.cc
 *
 *    Description:  primary particle position, direction, energy etc.
 *
 *        Version:  1.0
 *        Created:  11.10.2009 14:43:03
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <globals.hh>
#include <Randomize.hh>
#include "CexmcPrimaryGeneratorAction.hh"
#include "CexmcPrimaryGeneratorActionMessenger.hh"
#include "CexmcParticleGun.hh"
#include "CexmcCommon.hh"


CexmcPrimaryGeneratorAction::CexmcPrimaryGeneratorAction(
                                    CexmcPhysicsManager *  physicsManager ) :
    particleGun( NULL ), fwhmPosX( 0 ), fwhmPosY( 0 ), fwhmDirX( 0 ),
    fwhmDirY( 0 ), messenger( NULL )
{
    particleGun = new CexmcParticleGun( physicsManager );
    messenger = new CexmcPrimaryGeneratorActionMessenger( this );
}


CexmcPrimaryGeneratorAction::~CexmcPrimaryGeneratorAction()
{
    delete particleGun;
    delete messenger;
}


void  CexmcPrimaryGeneratorAction::GeneratePrimaries( G4Event *  event )
{
    particleGun->PrepareForNewEvent();

    const G4ThreeVector &  origPos( particleGun->GetOrigPosition() );
    const G4ThreeVector &  origDir( particleGun->GetOrigDirection() );
    G4double               origMomentumAmp( particleGun->GetOrigMomentumAmp() );

    G4double       randPosX( G4RandGauss::shoot( origPos.x(),
                                            fwhmPosX * CexmcFwhmToStddev ) );
    G4double       randPosY( G4RandGauss::shoot( origPos.y(),
                                            fwhmPosY * CexmcFwhmToStddev ) );
    G4ThreeVector  newPos( randPosX, randPosY, origPos.z() );

    G4double       randAngleX( G4RandGauss::shoot( origDir.x(),
                                            fwhmDirX * CexmcFwhmToStddev ) );
    G4double       randAngleY( G4RandGauss::shoot( origDir.y(),
                                            fwhmDirY * CexmcFwhmToStddev ) );
    G4ThreeVector  newAngle( randAngleX, randAngleY, origDir.z() );

    G4double       newMomentumAmp( G4RandGauss::shoot( origMomentumAmp,
                                            fwhmMomentumAmp * origMomentumAmp *
                                            CexmcFwhmToStddev ) );

    particleGun->SetParticlePosition( newPos );
    particleGun->SetParticleMomentumDirection( newAngle );
    particleGun->SetParticleMomentum( newMomentumAmp );

    particleGun->GeneratePrimaryVertex( event );
}


CexmcParticleGun *  CexmcPrimaryGeneratorAction::GetParticleGun( void )
{
    return particleGun;
}

