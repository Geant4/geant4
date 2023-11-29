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

// Manager class for Root read ntuples.
// It implements functions specific to Root read ntuples.
//
// Author: Ivana Hrivnacova, 09/04/2014 (ivana@ipno.in2p3.fr)

#ifndef G4RootRNtupleManager_h
#define G4RootRNtupleManager_h 1

#include "G4TRNtupleManager.hh"
#include "globals.hh"

#include "tools/rroot/ntuple"

#include <string_view>
#include <utility>

class G4RootRFileManager;
struct G4RootRNtupleDescription;

class G4RootRNtupleManager : public G4TRNtupleManager<tools::rroot::ntuple>
{
  friend class G4RootAnalysisReader;

  public:
    explicit G4RootRNtupleManager(const G4AnalysisManagerState& state);
    G4RootRNtupleManager() = delete;
    ~G4RootRNtupleManager() override = default;

  private:
    // Set methods
    void SetFileManager(std::shared_ptr<G4RootRFileManager> fileManager);

    // Methods from the base class
    G4int ReadNtupleImpl(const G4String& ntupleName, const G4String& fileName,
      const G4String& dirName, G4bool isUserFileName) final;
    G4bool GetTNtupleRow(G4TRNtupleDescription<tools::rroot::ntuple>* ntupleDescription) final;

    // Static data members
    static constexpr std::string_view fkClass { "G4RootPNtupleManager" };

    // Data members
    std::shared_ptr<G4RootRFileManager>  fFileManager { nullptr };
};

inline void
G4RootRNtupleManager::SetFileManager(std::shared_ptr<G4RootRFileManager> fileManager)
{
  fFileManager = std::move(fileManager);
}

#endif

