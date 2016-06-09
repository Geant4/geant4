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
 *       Filename:  CexmcPhysicsManager.hh
 *
 *    Description:  interface for external access to physics aspects
 *
 *        Version:  1.0
 *        Created:  27.10.2009 23:10:31
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_PHYSICS_MANAGER_HH
#define CEXMC_PHYSICS_MANAGER_HH

#include <G4Types.hh>
#include <G4ThreeVector.hh>

class  G4ParticleDefinition;
class  G4Track;
class  G4StepPoint;
class  CexmcProductionModel;
class  CexmcSetup;
class  CexmcPhysicsManagerMessenger;


class  CexmcPhysicsManager
{
    public:
        CexmcPhysicsManager();

        virtual ~CexmcPhysicsManager();

    public:
        virtual CexmcProductionModel *  GetProductionModel( void ) = 0;

        virtual G4bool  IsStudiedProcessAllowed( void ) const = 0;

        virtual void    ResampleTrackLengthInTarget( const G4Track *  track,
                                    const G4StepPoint *  stepPoint = NULL ) = 0;

        virtual void    SetupConstructionHook( const CexmcSetup *  setup ) = 0;

    public:
        G4bool    OnlyBeamParticleCanTriggerStudiedProcess( void ) const;

        void      IncrementNumberOfTriggeredStudiedInteractions( void );

        void      ResetNumberOfTriggeredStudiedInteractions( void );

        G4double  GetProposedMaxIL( void ) const;

        void      SetMaxIL( const G4ThreeVector &  direction );

        void      SetMaxILCorrection( G4double  value );

        void      SetProposedMaxIL( G4double  value );

    protected:
        virtual void    CalculateBasicMaxIL(
                                        const G4ThreeVector &  direction ) = 0;

    protected:
        G4double  basicMaxIL;

        G4double  maxILCorrection;

        G4double  proposedMaxIL;

        G4int     numberOfTriggeredStudiedInteractions;

        G4bool    onlyBeamParticleCanTriggerStudiedProcess;

    private:
        CexmcPhysicsManagerMessenger *  messenger;
};


inline G4bool  CexmcPhysicsManager::OnlyBeamParticleCanTriggerStudiedProcess(
                                                                    void ) const
{
    return onlyBeamParticleCanTriggerStudiedProcess;
}


inline void  CexmcPhysicsManager::IncrementNumberOfTriggeredStudiedInteractions(
                                                                        void )
{
    ++numberOfTriggeredStudiedInteractions;
}


inline void  CexmcPhysicsManager::ResetNumberOfTriggeredStudiedInteractions(
                                                                        void )
{
    numberOfTriggeredStudiedInteractions = 0;
}


inline G4double  CexmcPhysicsManager::GetProposedMaxIL( void ) const
{
    return proposedMaxIL;
}


inline void  CexmcPhysicsManager::SetMaxIL( const G4ThreeVector &  direction )
{
    CalculateBasicMaxIL( direction );
    proposedMaxIL = basicMaxIL + maxILCorrection;
}


inline void  CexmcPhysicsManager::SetMaxILCorrection( G4double  value )
{
    maxILCorrection = value;
    proposedMaxIL = basicMaxIL + maxILCorrection;
}


inline void  CexmcPhysicsManager::SetProposedMaxIL( G4double  value )
{
    proposedMaxIL = value;
}


#endif

