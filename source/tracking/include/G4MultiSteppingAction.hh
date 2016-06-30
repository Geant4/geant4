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
//
// $Id: G4MultiSteppingAction.hh 90212 2016-01-27 18:33:12Z adotti $
//
//---------------------------------------------------------------
//
// G4MultiSteppingAction.hh
//
//   Created on: Jan 17, 2016
//       Author: adotti
//
//
// class description:
//     This class extends G4UserSteppingAction and allows multiple
//     user-defined tracking actions to be used in the same job.
//     The class is a vector of user-defined tracking actions.
//     This class owns and manages the dependent user-actions.
// Usage:
//     There is no need to explicitly use this class as long as the
//     user actions are set via G4UserActionInitialization::SetUserAction
//     that can be called several times. Explicitly this is what is happening:
//     In user-defined action initialization:
//      G4MultiSteppingAction* action = new G4MultiSteppingAction;
//      action->push_back( G4UserSteppingActionUPtr( new MyUserSteppingAction );
//      [... add as many as needed ...]
//      SetUserAction( action );
// ---------------------------------------------------------------
//

#ifndef G4MULTISTEPPINGACTION_HH_
#define G4MULTISTEPPINGACTION_HH_

#include "G4UserSteppingAction.hh"
#include <vector>
#include <memory>

using G4UserSteppingActionUPtr=std::unique_ptr<G4UserSteppingAction>;
using G4UserSteppingActionVector=std::vector<G4UserSteppingActionUPtr>;

class G4MultiSteppingAction : public G4UserSteppingAction , public G4UserSteppingActionVector
{
public:
  G4MultiSteppingAction() = default;
  virtual ~G4MultiSteppingAction() override = default;
  virtual void UserSteppingAction(const G4Step*) override;
  virtual void SetSteppingManagerPointer(G4SteppingManager* pValue) override;
};

#endif /* SOURCE_TRACKING_INCLUDE_G4MULTISTEPPINGACTION_HH_ */
