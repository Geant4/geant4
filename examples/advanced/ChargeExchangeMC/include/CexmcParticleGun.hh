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
 *       Filename:  CexmcParticleGun.hh
 *
 *    Description:  particle gun
 *
 *        Version:  1.0
 *        Created:  15.12.2009 00:41:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PARTICLE_GUN_HH
#define CEXMC_PARTICLE_GUN_HH

#include <G4ParticleGun.hh>
#include <G4ThreeVector.hh>
#include "CexmcPhysicsManager.hh"
#include "CexmcException.hh"

class  CexmcParticleGunMessenger;


class  CexmcParticleGun : public G4ParticleGun
{
    public:
        explicit CexmcParticleGun( CexmcPhysicsManager *  physicsManager,
                                   G4int  nmbOfParticles = 1 );

        ~CexmcParticleGun();

    public:
        void  PrepareForNewEvent( void );

    public:
        const G4ThreeVector &  GetOrigPosition( void ) const;

        const G4ThreeVector &  GetOrigDirection( void ) const;

        G4double               GetOrigMomentumAmp( void ) const;

        void  SetOrigPosition( const G4ThreeVector &  position,
                               G4bool  fromMessenger = true );

        void  SetOrigDirection( const G4ThreeVector &  direction,
                                G4bool  fromMessenger = true );

        void  SetOrigMomentumAmp( G4double  momentumAmp,
                                  G4bool  fromMessenger = true );

        void  SetBeamParticle( G4ParticleDefinition *  particleDefinition,
                               G4bool  fromMessenger = true );

    private:
        CexmcPhysicsManager *        physicsManager;

        G4ThreeVector                origPos;

        G4ThreeVector                origDir;

        G4double                     origMomentumAmp;

    private:
        CexmcParticleGunMessenger *  messenger;
};


inline void  CexmcParticleGun::PrepareForNewEvent( void )
{
    /* this will prevent G4ParticleGun spam about kinetic energy redefinition */
    particle_energy = 0.0;
    particle_momentum = 0.0;
}


inline const G4ThreeVector &  CexmcParticleGun::GetOrigPosition( void ) const
{
    return origPos;
}


inline const G4ThreeVector &  CexmcParticleGun::GetOrigDirection( void ) const
{
    return origDir;
}


inline G4double  CexmcParticleGun::GetOrigMomentumAmp( void ) const
{
    return origMomentumAmp;
}


inline void  CexmcParticleGun::SetOrigPosition(
                    const G4ThreeVector &  position, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    origPos = position;
}


inline void  CexmcParticleGun::SetOrigDirection(
                    const G4ThreeVector &  direction, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    origDir = direction;

    physicsManager->SetMaxIL( direction );
}


inline void  CexmcParticleGun::SetOrigMomentumAmp( G4double  momentumAmp,
                                                   G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    origMomentumAmp = momentumAmp;
}


inline void  CexmcParticleGun::SetBeamParticle(
            G4ParticleDefinition *  particleDefinition, G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    SetParticleDefinition( particleDefinition );
}


#endif

