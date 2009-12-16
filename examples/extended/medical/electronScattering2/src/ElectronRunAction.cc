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

#include "ElectronRunAction.hh"
#include "ElectronRun.hh"
#include "G4Run.hh"
#include <assert.h>

ElectronRunAction::ElectronRunAction(G4String &outputFile) {
  outputFileSpec = outputFile;
}

ElectronRunAction::~ElectronRunAction() {}

void ElectronRunAction::BeginOfRunAction(const G4Run*) {}

G4Run*  ElectronRunAction::GenerateRun()
{
  return new ElectronRun("MyDetector");
}

void ElectronRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout <<"Number of Events Processed:" <<aRun->GetNumberOfEvent() <<G4endl;

  const ElectronRun* theRun = dynamic_cast<const ElectronRun*>(aRun);
  assert (0 != theRun);

  theRun->DumpData(outputFileSpec);
}
