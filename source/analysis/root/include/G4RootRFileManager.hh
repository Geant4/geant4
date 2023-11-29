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

// The manager for Root file input operations.

// Author: Ivana Hrivnacova, 10/09/2014  (ivana@ipno.in2p3.fr)

#ifndef G4RootRFileManager_h
#define G4RootRFileManager_h 1

#include "G4VRFileManager.hh"
#include "G4RootRFileDef.hh"
#include "globals.hh"

#include <map>
#include <string_view>

class G4RootRFileManager : public G4VRFileManager
{
  public:
    explicit G4RootRFileManager(const G4AnalysisManagerState& state);
    G4RootRFileManager() = delete;
    ~G4RootRFileManager() override;

    G4String GetFileType() const final { return "Root"; }

    // Methods from base class
    void CloseFiles() final {}

    // Methods to manipulate input files
    virtual G4bool OpenRFile(const G4String& fileName,
                             G4bool isPerThread);

    // Get methods
    G4RootRFile* GetRFile(const G4String& fileName, G4bool isPerThread) const;

  private:
    // Static data members
    static constexpr std::string_view fkClass { "G4RootRFileManager" };

    // data members
    std::map<G4String, G4RootRFile*> fRFiles;
};

#endif
