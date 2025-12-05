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

#ifndef RunAction_h
#define RunAction_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "G4UserRunAction.hh"

#include "globals.hh"

class G4Run;
class G4IAEAphspWriterStack;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class RunAction : public G4UserRunAction
{
public:

  RunAction() = default;
  virtual ~RunAction() override;

  G4Run* GenerateRun() override;
  // A derived G4Run is needed to store IAEAphsp particles during local run
  // and to dump info into the IAEAphsp files using Run::Merge()

  // void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

  // Modifiers and setters
  void SetIAEAphspWriterStack(const G4String& namePrefix);
  void AddZphsp(const G4double val);


private:

  // IAEAphsp stack object for the local run
  G4IAEAphspWriterStack* fIAEAphspWriterStack = nullptr;

};

#endif
