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

// Extension of Root ntuple manager for MPI.
// The class handles ntuples on the collecting MPI ranks.
// This class is temporarily provided with g4mpi,
// it will be integrated in Geant4 analysis category in future.
//
// Author: Ivana Hrivnacova, 21/11/2018 (ivana@ipno.in2p3.fr)

#ifndef G4RootMpiNtupleManager_h
#define G4RootMpiNtupleManager_h 1

#include "G4RootNtupleManager.hh"

namespace tools {
class impi;    
}

class G4RootMpiNtupleManager : public G4RootNtupleManager
{
  friend class G4RootMpiAnalysisManager;
  friend class G4RootMpiMainNtupleManager;

  public:
    G4RootMpiNtupleManager(const G4AnalysisManagerState& state, 
                           G4bool rowWise, G4bool rowMode,
                           tools::impi* impi, G4int mpiSize);
    virtual ~G4RootMpiNtupleManager();

    virtual void CreateNtuplesFromBooking() final;
    virtual G4bool Merge() final;

   private:
    // MPI
    G4bool Send(G4int id, tools::wroot::ntuple* ntuple);
    G4bool InitializeRanks();
    G4bool WaitBuffer();

    tools::impi*  fImpi;
    std::vector<G4int>  fSlaveRanks;
    G4int fMainRank;
};    

#endif


