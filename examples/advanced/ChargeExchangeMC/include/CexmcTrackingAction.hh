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

