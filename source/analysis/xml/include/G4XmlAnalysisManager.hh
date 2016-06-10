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
// $Id: G4XmlAnalysisManager.hh 74257 2013-10-02 14:24:55Z gcosmo $

// The main manager for Xml analysis.
// It delegates most of functions to the object specific managers. 

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4XmlAnalysisManager_h
#define G4XmlAnalysisManager_h 1

#include "G4VAnalysisManager.hh"
#include "globals.hh"

#include "tools/histo/h1d" 
#include "tools/histo/h2d" 
#include "tools/waxml/ntuple"


class G4XmlFileManager;
class G4H1ToolsManager;
class G4H2ToolsManager;
class G4XmlNtupleManager;

class G4XmlAnalysisManager : public G4VAnalysisManager
{
  public:
    G4XmlAnalysisManager(G4bool isMaster = true);
    ~G4XmlAnalysisManager();
    
    // static methods
    static G4XmlAnalysisManager* Instance();

    // Access methods
    tools::histo::h1d* GetH1(G4int id, G4bool warn = true,
                             G4bool onlyIfActive = true) const;
    tools::histo::h2d* GetH2(G4int id, G4bool warn = true,
                             G4bool onlyIfActive = true) const;
    tools::waxml::ntuple* GetNtuple() const;
    tools::waxml::ntuple* GetNtuple(G4int ntupleId) const;

  protected:
    // virtual methods from base class
    virtual G4bool OpenFileImpl(const G4String& fileName);
    virtual G4bool WriteImpl();
    virtual G4bool CloseFileImpl(); 

  private:
    // static data members
    static G4XmlAnalysisManager* fgMasterInstance;
    static G4ThreadLocal G4XmlAnalysisManager* fgInstance;

    // methods
    G4bool WriteH1();
    G4bool WriteH2();
    G4bool WriteNtuple();
    G4bool CloseNtupleFiles();
    G4bool Reset();

    // data members
    G4H1ToolsManager*    fH1Manager;
    G4H2ToolsManager*    fH2Manager;
    G4XmlNtupleManager*  fNtupleManager;
    G4XmlFileManager*    fFileManager;

};

#include "G4XmlAnalysisManager.icc"

#endif

