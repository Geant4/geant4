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

// The helper class for histogram and profiles messengers.
// It implements reusable commands in /analysis/[hn|pn] directory.
//
// Author: Ivana Hrivnacova, 05/05/2015  (ivana@ipno.in2p3.fr)

#ifndef G4AnalysisMessengerHelper_h
#define G4AnalysisMessengerHelper_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>
#include <string_view>

class G4VAnalysisManager;
class G4UIdirectory;
class G4UIcommand;

class G4AnalysisMessengerHelper
{
  public:
    // types
    struct BinData {
      BinData() {}
      G4int    fNbins { 0 };
      G4double fVmin { 0. };
      G4double fVmax { 0. };
      G4String fSunit;
      G4String fSfcn;
      G4String fSbinScheme;
    };
    // types
    struct ValueData {
      ValueData() {}
      G4double fVmin { 0. };
      G4double fVmax { 0. };
      G4String fSunit;
      G4String fSfcn;
    };

  public:
    // Make available utility Update method
    friend class G4HnMessenger;

 public:
    explicit G4AnalysisMessengerHelper(const G4String& hnType);
    G4AnalysisMessengerHelper() = delete;
    ~G4AnalysisMessengerHelper() = default;

    // Methods to create commands
    std::unique_ptr<G4UIdirectory>  CreateHnDirectory() const;

    std::unique_ptr<G4UIcommand>  CreateGetCommand(
                                        G4UImessenger* messenger) const;

    std::unique_ptr<G4UIcommand>  CreateSetTitleCommand(
                                        G4UImessenger* messenger) const;
    std::unique_ptr<G4UIcommand>  CreateSetBinsCommand(const G4String& axis,
                                        G4UImessenger* messenger) const;
    std::unique_ptr<G4UIcommand>  CreateSetValuesCommand(const G4String& axis,
                                        G4UImessenger* messenger) const;
    std::unique_ptr<G4UIcommand>  CreateSetAxisCommand(const G4String& axis,
                                        G4UImessenger* messenger) const;
    std::unique_ptr<G4UIcommand>  CreateSetAxisLogCommand(const G4String& axis,
                                        G4UImessenger* messenger) const;

    // Methods to read command paremeters
    void GetBinData(BinData& data, std::vector<G4String>& parameters,
                    G4int& counter) const;
    void GetValueData(ValueData& data, std::vector<G4String>& parameters,
                    G4int& counter) const;

    // Set methods
    void SetHnType(const G4String& hnType);

    // warnings
    void WarnAboutParameters(G4UIcommand* command, G4int nofParameters) const;
    void WarnAboutSetCommands() const;

  private:
    // Methods
    G4String Update(const G4String& str, const G4String& axis = "") const;

    // Static data members
    static constexpr std::string_view fkClass { "G4AnalysisMessengerHelper" };

    // Data members
    G4String  fHnType;
};

// inline functions

inline void G4AnalysisMessengerHelper::SetHnType(const G4String& hnType)
{ fHnType = hnType; }

#endif

