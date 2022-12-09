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

// The manager for Csv file input operations.

// Author: Ivana Hrivnacova, 21/10/2014  (ivana@ipno.in2p3.fr)

#ifndef G4CsvRFileManager_h
#define G4CsvRFileManager_h 1

#include "G4VRFileManager.hh"
#include "globals.hh"

#include <fstream>
#include <map>
#include <string_view>

class G4AnalysisManagerState;

class G4CsvRFileManager : public G4VRFileManager
{
  public:
    explicit G4CsvRFileManager(const G4AnalysisManagerState& state);
    G4CsvRFileManager() = delete;
    ~G4CsvRFileManager() override;

    G4String GetFileType() const final { return "csv"; }

    // Methods from base class
    void CloseFiles() final {}

    // Methods to manipulate input files
    virtual G4bool OpenRFile(const G4String& fileName);

    // Specific methods for files per objects
    std::ifstream* GetRFile(const G4String& fileName) const;

   private:
    // Static data members
    static constexpr std::string_view fkClass { "G4CsvRFileManager" };

    // data members
    std::map<G4String, std::ifstream*> fRFiles;
};

#endif
