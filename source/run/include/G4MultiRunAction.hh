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
// $Id: G4MultiRunAction.hh 90212 2016-01-27 18:33:12Z adotti $
//
//---------------------------------------------------------------
//
// G4MultiRunAction.hh
//
//   Created on: Jan 17, 2016
//       Author: adotti
//
//
// class description:
//     This class extends G4UserRunAction and allows multiple
//     user-defined tracking actions to be used in the same job.
//     The class is a vector of user-defined tracking actions.
//     This class owns and manages the dependent user-actions.
// Usage:
//     There is no need to explicitly use this class as long as the
//     user actions are set via G4UserActionInitialization::SetUserAction
//     that can be called several times. Explicitly this is what is happening:
//     In user-defined action initialization:
//      G4MultiRunAction* action = new G4MultiRunAction;
//      action->push_back( G4UserRunActionUPtr( new MyUserRunAction );
//      [... add as many as needed ...]
//      SetUserAction( action );
// ---------------------------------------------------------------
//

#ifndef G4MULTIRUNACTION_HH_
#define G4MULTIRUNACTION_HH_

#include "G4UserRunAction.hh"
#include <vector>
#include <memory>

using G4UserRunActionUPtr=std::unique_ptr<G4UserRunAction>;
using G4UserRunActionVector=std::vector<G4UserRunActionUPtr>;


class G4MultiRunAction : public G4UserRunAction , public G4UserRunActionVector
{
public:
  G4MultiRunAction() = default;
  virtual ~G4MultiRunAction() = default;
  virtual G4Run* GenerateRun() override;
  virtual void BeginOfRunAction(const G4Run* aRun) override;
  virtual void EndOfRunAction(const G4Run* aRun) override;
  virtual void SetMaster(G4bool val=true) override;

};

#endif /* SOURCE_RUN_INCLUDE_G4MULTIRUNACTION_HH_ */
