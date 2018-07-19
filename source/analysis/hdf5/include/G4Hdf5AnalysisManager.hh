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
// $Id$

// The main manager for Hdf5 analysis.
// It delegates most of functions to the object specific managers. 

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5AnalysisManager_h
#define G4Hdf5AnalysisManager_h 1

#include "G4ToolsAnalysisManager.hh"
#include "globals.hh"

#include "tools/hdf5/ntuple"

#include <memory>

class G4Hdf5FileManager;
class G4Hdf5NtupleManager;

class G4Hdf5AnalysisManager : public G4ToolsAnalysisManager
{
  public:
    explicit G4Hdf5AnalysisManager(G4bool isMaster = true);
    ~G4Hdf5AnalysisManager();
    
    // static methods
    static G4Hdf5AnalysisManager* Instance();
    static G4bool IsInstance();

    // Access methods
    tools::hdf5::ntuple* GetNtuple() const;
    tools::hdf5::ntuple* GetNtuple(G4int ntupleId) const;

    // Iterators
    std::vector<tools::hdf5::ntuple*>::iterator BeginNtuple();
    std::vector<tools::hdf5::ntuple*>::iterator EndNtuple();
    std::vector<tools::hdf5::ntuple*>::const_iterator BeginConstNtuple() const;
    std::vector<tools::hdf5::ntuple*>::const_iterator EndConstNtuple() const;
    
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
    static G4Hdf5AnalysisManager* fgMasterInstance;
    static G4ThreadLocal G4Hdf5AnalysisManager* fgInstance;

    // methods
    template <typename T>
    G4bool WriteHn(const std::vector<T*>& htVector,
                   const std::vector<G4HnInformation*>& hnVector,
                   const G4String& directoryName,
                   const G4String& hnType);
    template <typename T>
    G4bool WritePn(const std::vector<T*>& htVector,
                   const std::vector<G4HnInformation*>& hnVector,
                   const G4String& directoryName,
                   const G4String& hnType);
    G4bool WriteDirectory(); 
    G4bool WriteH1(); 
    G4bool WriteH2();
    G4bool WriteH3();
    G4bool WriteP1();
    G4bool WriteP2();
    G4bool Reset();

    // data members
    G4Hdf5NtupleManager*  fNtupleManager;
    std::shared_ptr<G4Hdf5FileManager>  fFileManager;
};

#include "G4Hdf5AnalysisManager.icc"

#endif

