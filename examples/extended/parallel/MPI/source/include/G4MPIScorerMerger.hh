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
#ifndef G4MPISCORERMERGER_HH
#define G4MPISCORERMERGER_HH
#include "G4ScoringManager.hh"
#include <vector>
#include <mpi.h>
#include "G4MPImanager.hh"

class G4MPIScorerMerger {
public:
  G4MPIScorerMerger( G4ScoringManager* mgr, 
                     G4int destination =  G4MPImanager::kRANK_MASTER,
                     G4int verbosity = 0 );
  virtual ~G4MPIScorerMerger() { clear(); }
  void SetDestinationRank( G4int i ) { destinationRank = i; }
  G4int GetDestinationRank() const { return destinationRank; }
  void SetScoringManager( G4ScoringManager* mgr ) { scoringManager = mgr; }
  G4ScoringManager* GetScoringManager() const { return scoringManager; }
  G4int GetCommSize() const { return commSize; }
  virtual void Merge();
  void SetVerbosity( G4int ver ) { verbose = ver; }
  G4int GetVerbosity() const { return verbose; }

protected:
  //Internal MPI-friendly format
  //for a single MeshScoreMap
  struct convMap_t {
    G4String name;
    G4int numElems;
    G4int* indexes;
    G4double* values;
  };
  virtual convMap_t* convertMap( const G4String& mapName ,
                                 G4THitsMap<double>* map ) const;
  virtual void convertMesh( const G4VScoringMesh* mesh );
  void clear();
  std::vector<convMap_t*> convertedMesh;
  G4int meshID;

  G4ScoringManager* scoringManager;
  G4int commSize;
  G4int destinationRank;
  MPI::Intracomm COMM_G4COMMAND_;
  G4int verbose;

  virtual void SendOneMesh();
  virtual void ReceiveOneMesh();
  virtual void MergeOneMesh();
  friend std::ostream& operator<<(std::ostream& os ,  const convMap_t& cnv );
};

#endif //G4MPISCORERMERGER_HH

