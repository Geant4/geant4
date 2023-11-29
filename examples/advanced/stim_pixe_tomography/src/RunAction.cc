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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4EmCalculator.hh"
#include "G4Run.hh"

#include "Run.hh"
#include "RunActionMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(G4bool isP)
  : G4UserRunAction(),
    fRun(nullptr),
    isPIXE(isP),
    fNbProjections(1),
    fNbSlices(1),
    fNbPixels(1),
    fInterFlag(false),
    fInterruptedProjectionIndex(0)

{
  fRunMessenger = new RunActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{
  fRun = new Run(isPIXE, fNbProjections, fNbSlices, fNbPixels, fInterruptedProjectionIndex);
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  if (IsMaster()) {
    if (fInterFlag) {
      G4cout << "####Resume simulation####" << G4endl;
    }

    G4cout << "Scan information: " << G4endl << "Projection_index " << fRun->GetCurrentProjection()
           << "; Slice_index " << fRun->GetCurrentSlice() << "; Pixel_index "
           << fRun->GetCurrentPixel() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  if (isMaster) fRun->EndOfRun();
}
void RunAction::SetScanParameters(G4int nbProjections, G4int nbSlices, G4int nbPixels)
{
  fNbProjections = nbProjections;
  fNbSlices = nbSlices;
  fNbPixels = nbPixels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetInterruptionFlag(G4bool inter)
{
  fInterFlag = inter;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetResumeProjectionIndex(G4int proj)
{
  if (fInterFlag) fInterruptedProjectionIndex = proj;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
