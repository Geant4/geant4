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
// Merge G4analysis histogram objects via MPI
//
// History:
// Jun 27, 2015 : Ivana Hrivnacova - new implementation using g4analysis

#include "G4MPIhistoMerger.hh"

#include "toolx/mpi/hmpi"

#include "G4MPImanager.hh"
#include "G4VAnalysisManager.hh"
#include "G4ios.hh"

#include <mpi.h>

G4MPIhistoMerger::G4MPIhistoMerger()
  : manager(0), destination(G4MPImanager::kRANK_MASTER), verboseLevel(0)
{}

G4MPIhistoMerger::G4MPIhistoMerger(G4VAnalysisManager* m, G4int dest, G4int v)
  : manager(m), destination(dest), verboseLevel(v)
{}

void G4MPIhistoMerger::Merge()
{
  if (verboseLevel > 0) {
    G4cout << "Starting merging of histograms" << G4endl;
  }

  const MPI_Comm* parentComm = G4MPImanager::GetManager()->GetComm();
  MPI_Comm comm;
  MPI_Comm_dup(*parentComm, &comm);

  G4bool verbose = (verboseLevel > 1);
  G4int tag = G4MPImanager::kTAG_HISTO;
  // const MPI::Intracomm* comm = &COMM_G4COMMAND_;
  toolx::mpi::hmpi* hmpi = new toolx::mpi::hmpi(G4cout, destination, tag, comm, verbose);
  if (!manager->Merge(hmpi)) {
    G4cout << " Merge FAILED" << G4endl;
  }

  delete hmpi;

  if (verboseLevel > 0) {
    G4cout << "End merging of histograms" << G4endl;
  }
  MPI_Comm_free(&comm);
}
