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


CexmcParticleGun::CexmcParticleGun( CexmcPhysicsManager *  physicsManager_,
                                    G4int  nmbOfParticles ) :
    G4ParticleGun( nmbOfParticles ), physicsManager( physicsManager_ ),
    origPos( 0, 0, 0 ), origDir( 0, 0, 0 ), origMomentumAmp( 0 ),
    messenger( NULL )
{
    messenger = new CexmcParticleGunMessenger( this );
}


CexmcParticleGun::~CexmcParticleGun()
{
    delete messenger;
}

