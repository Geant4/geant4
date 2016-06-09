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
 *       Filename:  CexmcTrackingAction.hh
 *
 *    Description:  tracking action
 *
 *        Version:  1.0
 *        Created:  22.11.2009 17:08:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACKING_ACTION_HH
#define CEXMC_TRACKING_ACTION_HH

#include <G4UserTrackingAction.hh>
#include "CexmcCommon.hh"

class  G4Track;
class  G4LogicalVolume;
class  G4ParticleDefinition;
class  CexmcPhysicsManager;


class  CexmcTrackingAction : public G4UserTrackingAction
{
    public:
        explicit CexmcTrackingAction( CexmcPhysicsManager *  physicsManager );

    public:
        void  PreUserTrackingAction( const G4Track *  track );

        void  BeginOfEventAction( void );

    private:
        void  ResetOutputParticleTrackId( void );

        void  ResetOutputParticleDecayProductCopyNumber( void );

        void  SetupIncidentParticleTrackInfo( const G4Track *  track );

    private:
        CexmcPhysicsManager *    physicsManager;

        const G4LogicalVolume *  targetVolume;

        G4int                    outputParticleTrackId;

        G4int                    outputParticleDecayProductCopyNumber;

    private:
        G4ParticleDefinition *   incidentParticle;

        G4ParticleDefinition *   outputParticle;

        G4ParticleDefinition *   nucleusOutputParticle;
};


inline void  CexmcTrackingAction::ResetOutputParticleTrackId( void )
{
    outputParticleTrackId = CexmcInvalidTrackId;
}


inline void  CexmcTrackingAction::ResetOutputParticleDecayProductCopyNumber(
                                                                        void )
{
    outputParticleDecayProductCopyNumber = 0;
}


inline void  CexmcTrackingAction::BeginOfEventAction( void )
{
    ResetOutputParticleTrackId();
    ResetOutputParticleDecayProductCopyNumber();
}


#endif

