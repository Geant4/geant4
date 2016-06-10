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
// $Id: G4CsvAnalysisManager.hh 74257 2013-10-02 14:24:55Z gcosmo $

// The main manager for Csv analysis.
// It delegates most of functions to the object specific managers. 

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4CsvAnalysisManager_h
#define G4CsvAnalysisManager_h 1


#include "G4VAnalysisManager.hh"
#include "G4CsvNtupleDescription.hh"
#include "globals.hh"

#include "tools/wcsv_ntuple"

#include <fstream>
#include <vector>
#include <map>

class G4H1DummyManager;
class G4H2DummyManager;
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
    tools::wcsv::ntuple* GetNtuple() const;
    tools::wcsv::ntuple* GetNtuple(G4int ntupleId) const;

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
    G4bool CloseNtupleFiles();
 
    // data members
    //
    G4H1DummyManager*   fH1Manager;
    G4H2DummyManager*   fH2Manager;
    G4CsvNtupleManager* fNtupleManager;
    G4CsvFileManager*   fFileManager;
};

#include "G4CsvAnalysisManager.icc"

#endif

