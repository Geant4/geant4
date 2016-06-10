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
// $Id: G4CsvAnalysisManager.hh 85025 2014-10-23 09:57:57Z gcosmo $

// The main manager for Csv analysis.
// It delegates most of functions to the object specific managers. 

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4CsvAnalysisManager_h
#define G4CsvAnalysisManager_h 1


#include "G4VAnalysisManager.hh"
#include "G4CsvNtupleDescription.hh"
#include "globals.hh"

#include "tools/histo/h1d" 
#include "tools/histo/h2d" 
#include "tools/histo/h3d" 
#include "tools/histo/p1d" 
#include "tools/histo/p2d" 
#include "tools/wcsv_ntuple"

#include <fstream>
#include <vector>
#include <map>

class G4H1ToolsManager;
class G4H2ToolsManager;
class G4H3ToolsManager;
class G4P1ToolsManager;
class G4P2ToolsManager;
class G4CsvFileManager;
class G4CsvNtupleManager;

class G4CsvAnalysisManager : public G4VAnalysisManager
{
  public:
    G4CsvAnalysisManager(G4bool isMaster = true);
    ~G4CsvAnalysisManager();
    
    // static methods
    static G4CsvAnalysisManager* Instance();

    // Access methods
    tools::histo::h1d* GetH1(G4int id, G4bool warn = true,
                             G4bool onlyIfActive = true) const;
    tools::histo::h2d* GetH2(G4int id, G4bool warn = true,
                             G4bool onlyIfActive = true) const;
    tools::histo::h3d*  GetH3(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;
    tools::histo::p1d*  GetP1(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;
    tools::histo::p2d*  GetP2(G4int id, G4bool warn = true,
                              G4bool onlyIfActive = true) const;

    tools::wcsv::ntuple* GetNtuple() const;
    tools::wcsv::ntuple* GetNtuple(G4int ntupleId) const;

    // Iterators
    std::vector<tools::histo::h1d*>::iterator BeginH1();
    std::vector<tools::histo::h1d*>::iterator EndH1();
    std::vector<tools::histo::h1d*>::const_iterator BeginConstH1() const;
    std::vector<tools::histo::h1d*>::const_iterator EndConstH1() const;
    
    std::vector<tools::histo::h2d*>::iterator BeginH2();
    std::vector<tools::histo::h2d*>::iterator EndH2();
    std::vector<tools::histo::h2d*>::const_iterator BeginConstH2() const;
    std::vector<tools::histo::h2d*>::const_iterator EndConstH2() const;
    
    std::vector<tools::histo::h3d*>::iterator BeginH3();
    std::vector<tools::histo::h3d*>::iterator EndH3();
    std::vector<tools::histo::h3d*>::const_iterator BeginConstH3() const;
    std::vector<tools::histo::h3d*>::const_iterator EndConstH3() const;
    
    std::vector<tools::histo::p1d*>::iterator BeginP1();
    std::vector<tools::histo::p1d*>::iterator EndP1();
    std::vector<tools::histo::p1d*>::const_iterator BeginConstP1() const;
    std::vector<tools::histo::p1d*>::const_iterator EndConstP1() const;
    
    std::vector<tools::histo::p2d*>::iterator BeginP2();
    std::vector<tools::histo::p2d*>::iterator EndP2();
    std::vector<tools::histo::p2d*>::const_iterator BeginConstP2() const;
    std::vector<tools::histo::p2d*>::const_iterator EndConstP2() const;

    std::vector<tools::wcsv::ntuple*>::iterator BeginNtuple();
    std::vector<tools::wcsv::ntuple*>::iterator EndNtuple();
    std::vector<tools::wcsv::ntuple*>::const_iterator BeginConstNtuple() const;
    std::vector<tools::wcsv::ntuple*>::const_iterator EndConstNtuple() const;
    
    // Csv format specific option
    void SetIsCommentedHeader(G4bool isCommentedHeader);
    void SetIsHippoHeader(G4bool isHippoHeader);

  protected:
    // virtual methods from base class
    virtual G4bool OpenFileImpl(const G4String& fileName);
    virtual G4bool WriteImpl();
    virtual G4bool CloseFileImpl(); 

  private:
    // static data members
    //
    static G4CsvAnalysisManager* fgMasterInstance;
    static G4ThreadLocal G4CsvAnalysisManager* fgInstance;
    
    // methods
    //
    G4bool WriteH1();
    G4bool WriteH2();
    G4bool WriteH3();
    G4bool WriteP1();
    G4bool WriteP2();
    G4bool CloseNtupleFiles();
    G4bool Reset();
 
    // data members
    //
    G4H1ToolsManager*   fH1Manager;
    G4H2ToolsManager*   fH2Manager;
    G4H3ToolsManager*   fH3Manager;
    G4P1ToolsManager*   fP1Manager;
    G4P2ToolsManager*   fP2Manager;
    G4CsvNtupleManager* fNtupleManager;
    G4CsvFileManager*   fFileManager;
};

#include "G4CsvAnalysisManager.icc"

#endif

