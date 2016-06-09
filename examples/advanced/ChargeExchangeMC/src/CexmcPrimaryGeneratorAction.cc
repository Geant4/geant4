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

