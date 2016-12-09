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
#include <memory>
#include <utility>
#include <mpi.h>
#include "G4MPImanager.hh"

//typedef G4THitsMap<G4double> HitMap;
typedef G4THitsMap<G4StatDouble> HitStatDoubleMap;

// This class allows for merging over MPI two command line scorers
// via MPI
class G4MPIscorerMerger {
public:
  G4MPIscorerMerger();
  G4MPIscorerMerger( G4ScoringManager* mgr, 
                     G4int destination =  G4MPImanager::kRANK_MASTER,
                     G4int verbosity = 0 );
  virtual ~G4MPIscorerMerger();

  //Get/set methods
  void SetDestinationRank( G4int i ) { destinationRank = i; }
  void SetScoringManager( G4ScoringManager* mgr ) { scoringManager = mgr; }
  void SetVerbosity( G4int ver ) { verbose = ver; }

  //Main Interface: call this method to merge all results to rank0
  void Merge();

protected:
  void SetupOutputBuffer(char* buff, G4int size, G4int position) {
    outputBuffer = buff;
    outputBufferSize=size;
    outputBufferPosition=position;
  }
  void DestroyBuffer() {
    delete[] outputBuffer;
    outputBuffer = nullptr;
    outputBufferSize=0;
    outputBufferPosition=0;
    ownsBuffer = false;
  }

  //! Pack all meshes into buffer
  void Pack(const G4ScoringManager*);
  void UnPackAndMerge(const G4ScoringManager*);

  //! Pack a single mesh
  void Pack(const G4VScoringMesh*);
  void UnPackAndMerge(G4VScoringMesh* );

  //! Pack a single score map
  //void Pack(const HitMap*);//Used When hits are <double>
  void Pack(const HitStatDoubleMap*);//Used when hits are statdouble
  //HitMap* UnPackHitMap(const G4String& detName, const G4String& colName);
  HitStatDoubleMap* UnPackHitStatDoubleMap(const G4String& detName, const G4String& colName);

  //Return size (in bytes) of the message needed to send the mesh
  G4int CalculatePackSize(const G4ScoringManager*) const;
  G4int CalculatePackSize(const G4VScoringMesh*) const;
  //G4int CalculatePackSize(const HitMap*) const;
  G4int CalculatePackSize(const HitStatDoubleMap* ) const;

protected:
  void Send(const unsigned int destination);
  void Receive(const unsigned int source);

private:
   char* outputBuffer;
   G4int outputBufferSize;
   G4int outputBufferPosition;
   long bytesSent;
   G4bool ownsBuffer;
  G4ScoringManager* scoringManager;
  unsigned int commSize;
  unsigned int destinationRank;
  MPI::Intracomm comm;
  G4int verbose;

};

#endif //G4MPISCORERMERGER_HH

