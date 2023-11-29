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

// Manager class for ntuple Xml file output.
//
// Author: Ivana Hrivnacova, 15/09/2020 (ivana@ipno.in2p3.fr)

#ifndef G4XmlNtupleFileManager_h
#define G4XmlNtupleFileManager_h 1

#include "G4VNtupleFileManager.hh"
#include "globals.hh"

#include "tools/waxml/ntuple"

#include <memory>
#include <utility>

class G4XmlFileManager;
class G4XmlNtupleManager;

class G4XmlNtupleFileManager : public G4VNtupleFileManager
{
  public:
    explicit G4XmlNtupleFileManager(const G4AnalysisManagerState& state);
    G4XmlNtupleFileManager() = delete;
    ~G4XmlNtupleFileManager() override = default;

    std::shared_ptr<G4VNtupleManager> CreateNtupleManager() override;

    // Methods to be performed at file management
    G4bool ActionAtOpenFile(const G4String& fileName) override;
    G4bool ActionAtWrite() override;
    G4bool ActionAtCloseFile() override;
    G4bool Reset() override;

    void SetFileManager(std::shared_ptr<G4XmlFileManager> fileManager);

    std::shared_ptr<G4XmlNtupleManager> GetNtupleManager() const;

  private:
    // Static data members
    static constexpr std::string_view fkClass { "G4XmlNtupleFileManager" };

    // Data members
    std::shared_ptr<G4XmlFileManager>  fFileManager { nullptr };
    std::shared_ptr<G4XmlNtupleManager>  fNtupleManager { nullptr };
};

// inline functions

inline void G4XmlNtupleFileManager::SetFileManager(
  std::shared_ptr<G4XmlFileManager> fileManager)
{
  fFileManager = std::move(fileManager);
}

inline std::shared_ptr<G4XmlNtupleManager> G4XmlNtupleFileManager::GetNtupleManager() const
{ return fNtupleManager; }

#endif

