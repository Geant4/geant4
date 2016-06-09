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
 *       Filename:  CexmcPrimaryGeneratorAction.hh
 *
 *    Description:  primary particle position, direction, energy etc.
 *
 *        Version:  1.0
 *        Created:  11.10.2009 14:54:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PRIMARY_GENERATOR_ACTION_HH
#define CEXMC_PRIMARY_GENERATOR_ACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>
#include "CexmcException.hh"

class  G4Event;
class  CexmcParticleGun;
class  CexmcPhysicsManager;
class  CexmcPrimaryGeneratorActionMessenger;


class  CexmcPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
    public:
        explicit CexmcPrimaryGeneratorAction(
                                        CexmcPhysicsManager *  physicsManager );

        ~CexmcPrimaryGeneratorAction();

    public:
        void      GeneratePrimaries( G4Event *  event );

    public:
        void      SetFwhmPosX( G4double  value, G4bool  fromMessenger = true );

        void      SetFwhmPosY( G4double  value, G4bool  fromMessenger = true );

        void      SetFwhmDirX( G4double  value, G4bool  fromMessenger = true );

        void      SetFwhmDirY( G4double  value, G4bool  fromMessenger = true );

        void      SetFwhmMomentumAmp( G4double  value,
                                      G4bool  fromMessenger = true );

        G4double  GetFwhmPosX( void ) const;

        G4double  GetFwhmPosY( void ) const;

        G4double  GetFwhmDirX( void ) const;

        G4double  GetFwhmDirY( void ) const;

        G4double  GetFwhmMomentumAmp( void ) const;

    public:
        CexmcParticleGun *  GetParticleGun( void );

    private:
        CexmcParticleGun *  particleGun;

        G4double            fwhmPosX;

        G4double            fwhmPosY;

        G4double            fwhmDirX;

        G4double            fwhmDirY;

        G4double            fwhmMomentumAmp;

    private:
        CexmcPrimaryGeneratorActionMessenger *  messenger;
};


inline void  CexmcPrimaryGeneratorAction::SetFwhmPosX( G4double  value,
                                                       G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    fwhmPosX = value;
}


inline void  CexmcPrimaryGeneratorAction::SetFwhmPosY( G4double  value,
                                                       G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    fwhmPosY = value;
}


inline void  CexmcPrimaryGeneratorAction::SetFwhmDirX( G4double  value,
                                                       G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    fwhmDirX = value;
}


inline void  CexmcPrimaryGeneratorAction::SetFwhmDirY( G4double  value,
                                                       G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    fwhmDirY = value;
}


inline void  CexmcPrimaryGeneratorAction::SetFwhmMomentumAmp( G4double  value,
                                                        G4bool  fromMessenger )
{
    if ( fromMessenger )
        ThrowExceptionIfProjectIsRead( CexmcCmdIsNotAllowed );

    fwhmMomentumAmp = value;
}


inline G4double  CexmcPrimaryGeneratorAction::GetFwhmPosX( void ) const
{
    return fwhmPosX;
}


inline G4double  CexmcPrimaryGeneratorAction::GetFwhmPosY( void ) const
{
    return fwhmPosY;
}


inline G4double  CexmcPrimaryGeneratorAction::GetFwhmDirX( void ) const
{
    return fwhmDirX;
}


inline G4double  CexmcPrimaryGeneratorAction::GetFwhmDirY( void ) const
{
    return fwhmDirY;
}


inline G4double  CexmcPrimaryGeneratorAction::GetFwhmMomentumAmp( void ) const
{
    return fwhmMomentumAmp;
}


#endif

