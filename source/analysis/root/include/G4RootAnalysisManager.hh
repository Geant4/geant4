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
class G4RootNtupleFileManager;

namespace tools {
namespace wroot {
class directory;    
}
}

class G4RootAnalysisManager : public  G4ToolsAnalysisManager
{
  friend class G4RootMpiAnalysisManager;

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
    virtual void SetNtupleMerging(G4bool mergeNtuples, 
                   G4int nofReducedNtupleFiles = 0) override;
    virtual void SetNtupleRowWise(G4bool rowWise, G4bool rowMode = true) override;
    virtual void SetBasketSize(unsigned int basketSize) override;
    virtual void SetBasketEntries(unsigned int basketEntries) override;

  protected:
    // virtual methods from base class
    virtual G4bool OpenFileImpl(const G4String& fileName) override;
    virtual G4bool WriteImpl() override;
    virtual G4bool CloseFileImpl(G4bool reset) override; 
    virtual G4bool IsOpenFileImpl() const final;

    // virtual functions (overriden in MPI implementation)
    virtual G4bool Reset();

  private:
    // static data members
    static G4RootAnalysisManager* fgMasterInstance;
    static G4ThreadLocal G4RootAnalysisManager* fgInstance;

    // methods
    template <typename T>
    G4bool WriteT(const std::vector<T*>& htVector,
                  const std::vector<G4HnInformation*>& hnVector,
                  const G4String& hnType);
    G4bool WriteH1();
    G4bool WriteH2();
    G4bool WriteH3();
    G4bool WriteP1();
    G4bool WriteP2();

    // data members 
    std::shared_ptr<G4RootFileManager> fFileManager;
    std::shared_ptr<G4RootNtupleFileManager> fNtupleFileManager;
};

#include "G4RootAnalysisManager.icc"

#endif
