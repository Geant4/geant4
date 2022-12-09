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

// Manager class for Xml ntuples
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4XmlNtupleManager_h
#define G4XmlNtupleManager_h 1

#include "G4TNtupleManager.hh"
#include "globals.hh"

#include "tools/waxml/ntuple"

#include <memory>
#include <string_view>
#include <utility>
#include <vector>

// Types alias
using XmlNtupleDescription = G4TNtupleDescription<tools::waxml::ntuple, std::ofstream>;

class G4XmlFileManager;

class G4XmlNtupleManager : public G4TNtupleManager<tools::waxml::ntuple,
                                                   std::ofstream>
{
  friend class G4XmlAnalysisManager;
  friend class G4XmlNtupleFileManager;

  public:
    explicit G4XmlNtupleManager(const G4AnalysisManagerState& state);
    G4XmlNtupleManager() = delete;
    ~G4XmlNtupleManager() override = default;

  private:
    // Functions specific to the output type
    //

    // Set methods
    void SetFileManager(std::shared_ptr<G4XmlFileManager> fileManager);

    // Access to ntuple vector (needed for Write())
    const std::vector<XmlNtupleDescription*>& GetNtupleDescriptionVector() const;

    // Methods from the templated base class
    //
    void CreateTNtupleFromBooking(XmlNtupleDescription* ntupleDescription) final;

    void FinishTNtuple(XmlNtupleDescription* ntupleDescription, G4bool fromBooking) final;

    // Static data members
    static constexpr std::string_view fkClass { "G4XmNtupleManager" };

    // Data members
    std::shared_ptr<G4XmlFileManager>  fFileManager { nullptr };
};

// inline functions

inline void
G4XmlNtupleManager::SetFileManager(std::shared_ptr<G4XmlFileManager> fileManager)
{
  fFileManager = std::move(fileManager);
}

inline const std::vector<G4TNtupleDescription<tools::waxml::ntuple, std::ofstream>*>&
G4XmlNtupleManager::GetNtupleDescriptionVector() const
{ return fNtupleDescriptionVector; }


#endif

