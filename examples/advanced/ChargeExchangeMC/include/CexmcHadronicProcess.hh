/*
 * =============================================================================
 *
 *       Filename:  CexmcHadronicProcess.hh
 *
 *    Description:  hadronic process with production model
 *
 *        Version:  1.0
 *        Created:  31.10.2009 23:44:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_HADRONIC_PROCESS_HH
#define CEXMC_HADRONIC_PROCESS_HH

#include <G4HadronicProcess.hh>
#include <G4Nucleus.hh>
#include "CexmcCommon.hh"

class  G4VParticleChange;
class  G4ParticleDefinition;
class  G4Track;
class  G4Step;
class  G4Material;
class  G4HadronicInteraction;
class  CexmcProductionModel;


class  CexmcHadronicProcess : public G4HadronicProcess
{
    public:
        explicit CexmcHadronicProcess(
                        const G4String &  name = CexmcStudiedProcessLastName );

        ~CexmcHadronicProcess();

    public:
        G4VParticleChange *  PostStepDoIt( const G4Track &  track,
                                           const G4Step &  step );

        G4bool  IsApplicable( const G4ParticleDefinition &  particle );

    public:
        void  RegisterProductionModel( CexmcProductionModel *  model );

    private:
        void  CalculateTargetNucleus( const G4Material *  material );

        void  FillTotalResult( G4HadFinalState *  hadFinalState,
                               const G4Track &  track );

    private:
        CexmcProductionModel *   productionModel;

        G4HadronicInteraction *  interaction;

    private:
        G4ParticleChange *      theTotalResult;

        G4Nucleus               targetNucleus;

        G4bool                  isInitialized;
};


#endif

