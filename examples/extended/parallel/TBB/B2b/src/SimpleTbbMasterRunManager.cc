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
//
// Class description:
//
//    This class implements the worker model run manager for TBB based
//    application.
//    It is instantiated by tbbUserWorkerInitialization and used by
//    SimpleTbbMasterRunManager.
//    See G4WorkerRunManager for documentation of methods relative to
//    base class. Only class specific methods are documented here.
//
// Equivalent in traditional MT:
//    G4WorkerRunManager
//
// History:
//    Oct 31st, 2013 A. Dotti - First Implementation

#include "SimpleTbbMasterRunManager.hh"
#include "tbbTask.hh"

SimpleTbbMasterRunManager::SimpleTbbMasterRunManager() : 
    G4MTRunManager(),
    theTasks(static_cast<tbb::task_list*>(nullptr)),
    nEvtsPerTask(1)
{
}

SimpleTbbMasterRunManager::~SimpleTbbMasterRunManager()
{
}

void SimpleTbbMasterRunManager::TerminateWorkers()
{
    //For TBB based example this should be empty
}

void SimpleTbbMasterRunManager::CreateAndStartWorkers()
{
    //Instead of pthread based workers, create tbbTask
    G4int ntasks = numberOfEventToBeProcessed/nEvtsPerTask;
    G4int remn = numberOfEventToBeProcessed % nEvtsPerTask;
    for ( G4int nt = 0 ; nt < ntasks ; ++nt )
    {
        G4int evts= nEvtsPerTask;
        if ( nt == ntasks - 1 ) evts+=remn;
        CreateTask(nt,evts);
    }
}

void SimpleTbbMasterRunManager::CreateTask(G4int id,G4int evts)
{
    tbbTask& task = * new(tbb::task::allocate_root())
    tbbTask( id , nullptr , evts ); //Add output for merging
    theTasks->push_back( task );
}

void SimpleTbbMasterRunManager::RunTermination()
{
    // Reduce results ....
    G4MTRunManager::RunTermination();
}

