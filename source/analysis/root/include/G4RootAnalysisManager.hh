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
// $Id: G4RootAnalysisManager.hh 106985 2017-10-31 10:07:18Z gcosmo $

// The main manager for Root analysis.
// It delegates most of functions to the object specific managers. 

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4RootAnalysisManager_h
#define G4RootAnalysisManager_h 1

#include "G4ToolsAnalysisManager.hh"
#include "globals.hh"

#include "tools/wroot/ntuple"
#include "tools/histo/hmpi"

#include <memory>

class G4RootFileManager;
class G4RootNtupleManager;
class G4RootMainNtupleManager;
class G4RootPNtupleManager;

namespace tools {
namespace wroot {
class directory;    
}
}

enum class G4NtupleMergeMode {
  kNone,
  kMain,
  kSlave
};

class G4RootAnalysisManager : public  G4ToolsAnalysisManager
{
  public:
    explicit G4RootAnalysisManager(G4bool isMaster = true);
    virtual ~G4RootAnalysisManager();
    
    // static methods
    static G4RootAnalysisManager* Instance();
    static G4bool IsInstance();

    // Access methods
    tools::wroot::ntuple* GetNtuple() const;
    tools::wroot::ntuple* GetNtuple(G4int ntupleId) const;

    // Iterators
    std::vector<tools::wroot::ntuple*>::iterator BeginNtuple();
    std::vector<tools::wroot::ntuple*>::iterator EndNtuple();
    std::vector<tools::wroot::ntuple*>::const_iterator BeginConstNtuple() const;
    std::vector<tools::wroot::ntuple*>::const_iterator EndConstNtuple() const;

    // MT/MPI
    void SetNtupleMerging(G4bool mergeNtuples, 
                          G4int nofReducedNtupleFiles = 0,
                          G4bool rowWise = 0,
                          unsigned int basketSize = fgkDefaultBasketSize);

  protected:
    // virtual methods from base class
    virtual G4bool OpenFileImpl(const G4String& fileName) final;
    virtual G4bool WriteImpl() final;
    virtual G4bool CloseFileImpl() final; 
    virtual G4bool IsOpenFileImpl() const final;

  private:
    // constants
    static constexpr unsigned int fgkDefaultBasketSize = 32000;

    // static data members
    static G4RootAnalysisManager* fgMasterInstance;
    static G4ThreadLocal G4RootAnalysisManager* fgInstance;    

    // methods
    void SetNtupleMergingMode(G4bool mergeNtuples, G4int nofNtupleFiles);
    void ClearNtupleManagers();
    void CreateNtupleManagers();
    G4int  GetNtupleFileNumber();

    template <typename T>
    G4bool WriteT(const std::vector<T*>& htVector,
                  const std::vector<G4HnInformation*>& hnVector,
                  tools::wroot::directory* directory,
                  const G4String& hnType);
    G4bool WriteH1();
    G4bool WriteH2();
    G4bool WriteH3();
    G4bool WriteP1();
    G4bool WriteP2();
    G4bool WriteNtuple();
    G4bool Reset();

    // data members 
    G4int   fNofNtupleFiles;
    G4bool  fNtupleRowWise;
    G4NtupleMergeMode      fNtupleMergeMode;
    G4RootNtupleManager*   fNtupleManager; 
    G4RootPNtupleManager*  fSlaveNtupleManager;
    std::shared_ptr<G4RootFileManager> fFileManager;
};

#include "G4RootAnalysisManager.icc"

#endif
