//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
