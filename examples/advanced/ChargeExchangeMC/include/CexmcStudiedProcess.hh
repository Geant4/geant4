/*
 * =============================================================================
 *
 *       Filename:  CexmcStudiedProcess.hh
 *
 *    Description:  studied process in the target
 *
 *        Version:  1.0
 *        Created:  26.10.2009 20:41:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_STUDIED_PROCESS_HH
#define CEXMC_STUDIED_PROCESS_HH

#include <G4WrapperProcess.hh>
#include <G4ProcessType.hh>

class  G4VParticleChange;
class  CexmcPhysicsManager;


class  CexmcStudiedProcess : public G4WrapperProcess
{
    public:
        explicit  CexmcStudiedProcess( CexmcPhysicsManager *  physicsManager,
                                   G4ProcessType  processType = fUserDefined );

    public:
        G4double  PostStepGetPhysicalInteractionLength( const G4Track &  track,
                    G4double  previousStepSize, G4ForceCondition *  condition );

        G4VParticleChange *  PostStepDoIt( const G4Track &  track,
                                           const G4Step &  step );

    private:
        CexmcPhysicsManager *  physicsManager;
};


#endif

