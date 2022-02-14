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
#ifndef RUNACTION_HH
#define RUNACTION_HH

#include "G4UserRunAction.hh"
#include "G4String.hh"

class G4Run;
class EventAction;
class G4GenericMessenger;

/**
 * @brief Run action
 *
 * Creates the output file with ntuple.
 * Output name can be changed with UI command "/HGCalTestbeam/output/file
 *  <NAME>"
 *
 */

class RunAction : public G4UserRunAction {
public:
  explicit RunAction(EventAction *);
  virtual ~RunAction();

  virtual void BeginOfRunAction(const G4Run *);
  virtual void EndOfRunAction(const G4Run *);

private:
  /// Pointer to the event action to retrieve vectors
  EventAction *fEventAction;
  /// Name of the output file
  G4String fOutputFileDir;
  /// Pointer to the command messenger
  G4GenericMessenger *fMessenger;
};

#endif /* RUNACTION_HH */
