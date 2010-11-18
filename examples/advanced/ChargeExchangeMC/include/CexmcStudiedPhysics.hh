/*
 * =============================================================================
 *
 *       Filename:  CexmcStudiedPhysics.hh
 *
 *    Description:  studied physics in the target
 *
 *        Version:  1.0
 *        Created:  18.10.2009 16:10:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_STUDIED_PHYSICS_HH
#define CEXMC_STUDIED_PHYSICS_HH

#include <G4VPhysicsConstructor.hh>
#include <G4ProcessManager.hh>
#include <G4ParticleDefinition.hh>
#include "CexmcStudiedProcess.hh"
#include "CexmcPhysicsManager.hh"
#include "CexmcProductionModel.hh"

class  G4VProcess;


template  < typename  Process >
class  CexmcStudiedPhysics : public G4VPhysicsConstructor
{
    public:
        explicit CexmcStudiedPhysics( CexmcPhysicsManager *  physicsManager );

        virtual ~CexmcStudiedPhysics();

    public:
        void                    ConstructParticle( void );

        void                    ConstructProcess( void );

    public:
        CexmcProductionModel *  GetProductionModel( void );

    protected:
        virtual void  ApplyInteractionModel( G4VProcess *  process );

    protected:
        CexmcPhysicsManager *   physicsManager;

        CexmcProductionModel *  productionModel;

    private:
        G4bool                  wasActivated;
};


template  < typename  Process >
CexmcStudiedPhysics< Process >::CexmcStudiedPhysics(
                                    CexmcPhysicsManager *  physicsManager ) :
    G4VPhysicsConstructor( "studiedPhysics" ), physicsManager( physicsManager ),
    productionModel( NULL ), wasActivated( false )
{
}


template  < typename  Process >
CexmcStudiedPhysics< Process >::~CexmcStudiedPhysics()
{
}


template  < typename  Process >
void  CexmcStudiedPhysics< Process >::ConstructParticle( void )
{
    if ( productionModel )
        productionModel->GetIncidentParticle();
}


template  <typename  Process >
void  CexmcStudiedPhysics< Process >::ConstructProcess( void )
{
    if ( wasActivated )
        return;

    wasActivated = true;

    Process *  process( new Process );

    CexmcStudiedProcess *  studiedProcess( new CexmcStudiedProcess(
                                physicsManager, process->GetProcessType() ) );

    ApplyInteractionModel( process );

    studiedProcess->RegisterProcess( process );

    G4ParticleDefinition *  particle( NULL );

    if ( productionModel )
        particle = productionModel->GetIncidentParticle();

    if ( particle )
    {
        G4ProcessManager *  processManager( particle->GetProcessManager() );
        processManager->AddDiscreteProcess( studiedProcess );
    }
}


template  < typename  Process >
CexmcProductionModel *  CexmcStudiedPhysics< Process >::GetProductionModel(
                                                                        void )
{
    return productionModel;
}


template  < typename  Process >
void  CexmcStudiedPhysics< Process >::ApplyInteractionModel( G4VProcess * )
{
}


#endif

