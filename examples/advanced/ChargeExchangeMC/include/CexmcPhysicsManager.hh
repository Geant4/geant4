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

