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
// $Id: G4XmlAnalysisManager.hh 92972 2015-09-23 14:36:03Z gcosmo $

// The main manager for Xml analysis.
// It delegates most of functions to the object specific managers. 

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4XmlAnalysisManager_h
#define G4XmlAnalysisManager_h 1

#include "G4ToolsAnalysisManager.hh"
#include "globals.hh"

#include "tools/waxml/ntuple"

#include <memory>

class G4XmlFileManager;
class G4XmlNtupleManager;

class G4XmlAnalysisManager : public G4ToolsAnalysisManager
{
  public:
    explicit G4XmlAnalysisManager(G4bool isMaster = true);
    ~G4XmlAnalysisManager();
    
    // static methods
    static G4XmlAnalysisManager* Instance();
    static G4bool IsInstance();

    // Access methods
    tools::waxml::ntuple* GetNtuple() const;
    tools::waxml::ntuple* GetNtuple(G4int ntupleId) const;

    // Iterators
    std::vector<tools::waxml::ntuple*>::iterator BeginNtuple();
    std::vector<tools::waxml::ntuple*>::iterator EndNtuple();
    std::vector<tools::waxml::ntuple*>::const_iterator BeginConstNtuple() const;
    std::vector<tools::waxml::ntuple*>::const_iterator EndConstNtuple() const;
    
  protected:
    // virtual methods from base class
    virtual G4bool OpenFileImpl(const G4String& fileName) final;
    virtual G4bool WriteImpl() final;
    virtual G4bool CloseFileImpl() final; 
    virtual G4bool IsOpenFileImpl() const final;

  private:
    // static data members
    static G4XmlAnalysisManager* fgMasterInstance;
    static G4ThreadLocal G4XmlAnalysisManager* fgInstance;

    // methods
    template <typename T>
    G4bool WriteT(const std::vector<T*>& htVector,
                  const std::vector<G4HnInformation*>& hnVector,
                  const G4String& directoryName,
                  const G4String& hnType);
    G4bool WriteH1();
    G4bool WriteH2();
    G4bool WriteH3();
    G4bool WriteP1();
    G4bool WriteP2();
    G4bool WriteNtuple();
    G4bool CloseNtupleFiles();
    G4bool Reset();

    // data members
    G4XmlNtupleManager*  fNtupleManager;
    std::shared_ptr<G4XmlFileManager>  fFileManager;
};

#include "G4XmlAnalysisManager.icc"

#endif

