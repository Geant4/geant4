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
// G4MultiEventAction
//
// Class description:
//
// This class extends G4UserEventAction and allows multiple user-defined
// tracking actions to be used in the same job.
// The class is a vector of user-defined tracking actions.
// This class owns and manages the dependent user-actions.
// Usage:
//   There is no need to explicitly use this class as long as the
//   user actions are set via G4UserActionInitialization::SetUserAction()
//   that can be called several times. Explicitly, this is what is happening:
//     In user-defined action initialization:
//     G4MultiEventAction* action = new G4MultiEventAction;
//     action->push_back( G4UserEventActionUPtr( new MyUserEventAction );
//     [... add as many as needed ...]
//     SetUserAction( action );

// Author: Andrea Dotti, SLAC - 17.01.2016
// --------------------------------------------------------------------
#ifndef G4MULTIEVENTACTION_HH
#define G4MULTIEVENTACTION_HH

#include "G4UserEventAction.hh"
#include <vector>
#include <memory>

using G4UserEventActionUPtr=std::unique_ptr<G4UserEventAction>;
using G4UserEventActionVector=std::vector<G4UserEventActionUPtr>;

class G4MultiEventAction : public G4UserEventAction
                         , public G4UserEventActionVector
{
  public:

    G4MultiEventAction() = default;
    ~G4MultiEventAction() override = default;
    void SetEventManager(G4EventManager* ) override;
    void BeginOfEventAction(const G4Event* ) override;
    void EndOfEventAction(const G4Event* ) override;
};

#endif
