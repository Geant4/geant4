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

#ifndef G4MPIHISTOMERGER_HH
#define G4MPIHISTOMERGER_HH

#include "G4MPImanager.hh"

class G4VAnalysisManager;

class G4MPIhistoMerger {
public:
  G4MPIhistoMerger();
  G4MPIhistoMerger(G4VAnalysisManager* mgr,
              G4int destination = G4MPImanager::kRANK_MASTER,
              G4int verbosity = 0);

  //Get/set methods
  void SetDestinationRank( G4int i ) { destination = i; }
  void SetScoringManager( G4VAnalysisManager* mgr ) { manager = mgr; }
  void SetVerbosity( G4int ver ) { verboseLevel = ver; }

  void Merge();

private:
  G4VAnalysisManager* manager;
  G4int destination;
  G4int verboseLevel;
};

#endif //G4MPIHISTOMERGERNEW_HH
