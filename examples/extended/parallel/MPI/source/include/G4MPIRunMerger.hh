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
#ifndef G4MPIRUNMERGER_HH
#define G4MPIRUNMERGER_HH
#include "G4Run.hh"
#include <mpi.h>
#include "G4MPImanager.hh"

class G4MPIRunMerger {
public:
  G4MPIRunMerger( const G4Run* aRun , 
                  G4int destination = G4MPImanager::kRANK_MASTER ,
                  G4int verbosity = 0);
  virtual ~G4MPIRunMerger() {}
  void SetRun( G4Run* r ) { run = r; }
  const G4Run* GetRun() const { return run; }
  void SetDestinationRank( G4int i ) { destinationRank = i; }
  G4int GetDestinationRank() const { return destinationRank; }
  G4int GetCommSize() const { return commSize; }
  virtual void Merge();
  void SetVerbosity( G4int ver ) { verbose = ver; }
  G4int GetVerbosity() const { return verbose; }
protected:
  virtual void Send();
  virtual void Receive(G4int rank);

  void SendDouble( G4double* val , G4int size=1);
  void SendInt( G4int* val , G4int size=1);
  void ReceiveDouble( G4int rank , G4double* val , G4int size=1);
  void ReceiveInt(G4int rank , G4int* val, G4int size=1);
  G4int destinationRank;
  G4Run* run;
  G4int commSize;
  MPI::Intracomm COMM_G4COMMAND_;
  G4int verbose;

};

#endif //G4MPIRUNMERGER_HH

