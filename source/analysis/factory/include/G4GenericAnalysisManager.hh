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

#ifndef G4GenericAnalysisManager_h
#define G4GenericAnalysisManager_h 1

#include "G4ToolsAnalysisManager.hh"
#include "G4THnManager.hh"
#include "G4Threading.hh"
#include "globals.hh"

#include <memory>

class G4GenericFileManager;
class G4HnInformation;
class G4VNtupleFileManager;

class G4GenericAnalysisManager : public  G4ToolsAnalysisManager
{
  friend class G4RootMpiAnalysisManager;

  public:
    explicit G4GenericAnalysisManager(G4bool isMaster = true);
    virtual ~G4GenericAnalysisManager();
    
    // static methods
    static G4GenericAnalysisManager* Instance();
    static G4bool IsInstance();

    // MT/MPI
    virtual void SetNtupleMerging(G4bool mergeNtuples,
                   G4int nofReducedNtupleFiles = 0) override;
    virtual void SetNtupleRowWise(G4bool rowWise, G4bool rowMode = true) override;
    virtual void SetBasketSize(unsigned int basketSize)  override;
    virtual void SetBasketEntries(unsigned int basketEntries)  override;

    // new functions
    G4bool Merge();
    
    // write in an extra file 
    G4bool WriteH1(G4int id, const G4String& fileName);
    G4bool WriteH2(G4int id, const G4String& fileName);
    G4bool WriteH3(G4int id, const G4String& fileName);
    G4bool WriteP1(G4int id, const G4String& fileName);
    G4bool WriteP2(G4int id, const G4String& fileName);

    // set default output type (backward compatibility)
    // this type will be used for file names without extension
    void SetDefaultFileType(const G4String& value);
    G4String GetDefaultFileType() const;

  protected:
    // virtual methods from base class
    virtual G4bool OpenFileImpl(const G4String& fileName) override;
    virtual G4bool WriteImpl() final;
    virtual G4bool CloseFileImpl(G4bool reset = true) override; 
    virtual G4bool IsOpenFileImpl() const final;
    virtual G4bool Reset() final;

  private:
    // methods
    void CreateNtupleFileManager(const G4String& fileName);

    // static data members
    static G4GenericAnalysisManager* fgMasterInstance;
    static G4ThreadLocal G4GenericAnalysisManager* fgInstance;

    // data members 
    std::shared_ptr<G4GenericFileManager> fFileManager;
    // add G4GenericNtupleManager
    // this class will be analogic to file managers but with ntuples
    std::shared_ptr<G4VNtupleFileManager> fNtupleFileManager;

    // data members
    G4bool  fIsNtupleMergingSet = false;
    G4int   fNofNtupleFiles = 0;
    G4bool  fMergeNtuples = false;
    G4bool  fNtupleRowWise = false;
    G4bool  fNtupleRowMode = true;
    unsigned int fBasketSize = fgkDefaultBasketSize;
    unsigned int fBasketEntries = fgkDefaultBasketEntries;
};

#include "G4GenericAnalysisManager.icc"

#endif
