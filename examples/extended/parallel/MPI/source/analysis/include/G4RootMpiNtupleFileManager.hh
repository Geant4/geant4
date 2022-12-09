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

// The Root analysis manager extension for MPI
// This class is temporarily provided with g4mpi,
// it will be integrated in Geant4 analysis category in future.
//
// Author: Ivana Hrivnacova, 21/11/2018 (ivana@ipno.in2p3.fr)

#ifndef G4RootMpiNtupleFileManager_h
#define G4RootMpiNtupleFileManager_h 1

#include "G4RootNtupleFileManager.hh"

class G4RootMpiNtupleManager;
class G4RootMpiPNtupleManager;

namespace tools {
class impi;  
}  

class G4RootMpiNtupleFileManager : public  G4RootNtupleFileManager
{
  public:
    explicit G4RootMpiNtupleFileManager(const G4AnalysisManagerState& state);
    virtual ~G4RootMpiNtupleFileManager();
    
    // MPI
    void SetMpiNtupleMerging(tools::impi* impi, 
                             G4int mpiRank, G4int mpiSize,
                             G4int nofReducedNtupleFiles = 0);

    // virtual methods from base class
    virtual G4bool ActionAtOpenFile(const G4String& fileName) final;
    virtual G4bool ActionAtWrite() final;
    virtual G4bool ActionAtCloseFile() final; 
    virtual G4bool Reset() final;

    virtual std::shared_ptr<G4VNtupleManager> CreateNtupleManager() override;

    std::shared_ptr<G4RootMpiNtupleManager> GetMpiNtupleManager() const;

  private:
    // methods
    void SetMpiNtupleMergingMode(G4int nofNtupleFiles);
    void CreateMpiNtupleManagers(tools::impi* impi, G4int mpiRank, G4int mpiSize);

    // data members 
    tools::impi* fImpi;
    G4int fMpiRank;
    G4int fMpiSize;
    // std::shared_ptr<G4RootMpiNtupleManager>  fMpiNtupleManager; 
    std::shared_ptr<G4RootMpiPNtupleManager>  fMpiSlaveNtupleManager;
    G4bool fNtupleBooked;
};

#endif
