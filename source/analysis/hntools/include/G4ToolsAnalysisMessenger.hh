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

// The messenger class for tools objects getters.
//
// Author: Ivana Hrivnacova, 29/10/2021  (ivana@ipno.in2p3.fr)

#ifndef G4ToolsAnalysisMessenger_h
#define G4ToolsAnalysisMessenger_h 1

#include "G4UImessenger.hh"
#include "G4THnManager.hh"
#include "globals.hh"

#include <memory>

class G4ToolsAnalysisManager;
class G4UIcommand;

class G4ToolsAnalysisMessenger : public G4UImessenger
{
  public:
    explicit G4ToolsAnalysisMessenger(G4ToolsAnalysisManager* manager);
    virtual ~G4ToolsAnalysisMessenger();

    // Methods
    G4String GetCurrentValue (G4UIcommand* command);
    virtual void SetNewValue(G4UIcommand* command, G4String value) final;

  private:
    // Methods
    template <typename HT>
    G4String GetHnAddress(G4int id, G4THnManager<HT>* htManager) const;

    // Data members
    G4ToolsAnalysisManager*  fManager; ///< Associated class
    std::unique_ptr<G4UIcommand>  fGetH1Cmd;
    std::unique_ptr<G4UIcommand>  fGetH2Cmd;
    std::unique_ptr<G4UIcommand>  fGetH3Cmd;
    std::unique_ptr<G4UIcommand>  fGetP1Cmd;
    std::unique_ptr<G4UIcommand>  fGetP2Cmd;
    G4String fH1Value;
    G4String fH2Value;
    G4String fH3Value;
    G4String fP1Value;
    G4String fP2Value;
};

// inline functions

template <typename HT>
inline
G4String G4ToolsAnalysisMessenger::GetHnAddress(
  G4int id, G4THnManager<HT>* htManager) const
{
  auto ht = htManager->GetT(id);
  if ( ht != nullptr ) {
    std::ostringstream os;
    os << static_cast<void*>(ht);
    return os.str();
  }
  return G4String();
}

#endif
