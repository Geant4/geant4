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
// Class for configuring G4analysis for merging ntuples via MPI
//
// Author: Ivana Hrivnacova, 21/11/2018 (ivana@ipno.in2p3.fr)

#include "G4MPIntupleMerger.hh"
#include "G4RootMpiAnalysisManager.hh"
#include "G4ios.hh"

#include "toolx/mpi/wrmpi"

#include <mpi.h>


G4MPIntupleMerger::G4MPIntupleMerger(G4int nofReducedNtupleFiles,
                                     G4bool rowWise, G4bool rowMode)
{
  // Configure MPI using  G4MPImanager

  // G4cout << "Start configure ntuple MPI merging" << G4endl;

  // Create MPI Root analysis manager
  G4bool isMaster = true;
  auto analysisManager
    = new G4RootMpiAnalysisManager(isMaster);
  analysisManager->SetVerboseLevel(1);
  // G4cout << "Start configure ntuple MPI merging" << G4endl;

  // Get communicator
  G4MPImanager* mpiManager = G4MPImanager::GetManager();
  G4int mpiRank = mpiManager->GetRank();
  G4int mpiSize = mpiManager->GetActiveSize();
  auto comm  = mpiManager->GetAllComm();

  // G4int tag = G4MPImanager::kTAG_NTUPLE;
  fWrmpi = new toolx::mpi::wrmpi(G4cout, *comm);

  analysisManager->SetMpiNtupleMerging(
    fWrmpi, mpiRank, mpiSize, nofReducedNtupleFiles);
  analysisManager->SetNtupleRowWise(rowWise, rowMode);

  // G4cout << "End configure ntuple MPI merging" << G4endl;
}

G4MPIntupleMerger::~G4MPIntupleMerger()
{
  delete fWrmpi;
}
