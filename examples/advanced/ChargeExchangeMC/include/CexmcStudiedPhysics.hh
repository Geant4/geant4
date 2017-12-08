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
#include "CexmcFakeCrossSectionData.hh"
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
                                    CexmcPhysicsManager *  physicsManager_ ) :
    G4VPhysicsConstructor( "studiedPhysics" ),
    physicsManager( physicsManager_ ), productionModel( NULL ),
    wasActivated( false )
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


template  < typename  Process >
void  CexmcStudiedPhysics< Process >::ConstructProcess( void )
{
    if ( wasActivated )
        return;

    wasActivated = true;

    Process *  process( new Process );

    process->AddDataSet( new CexmcFakeCrossSectionData );

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

