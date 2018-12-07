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
/// \file electromagnetic/TestEm3/src/EmAcceptance.cc
/// \brief Implementation of the Emeptance class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EmAcceptance.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmAcceptance::EmAcceptance()
 : fIsAccepted(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmAcceptance::~EmAcceptance()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmAcceptance::BeginOfAcceptance(const G4String& title, G4int stat)
{
  G4cout << "\n<<<<<ACCEPTANCE>>>>> " << stat << " events for " << title 
         << G4endl;
  fIsAccepted = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmAcceptance::EndOfAcceptance()
{
  G4String resume = "IS ACCEPTED";
  if (!fIsAccepted) resume = "IS NOT ACCEPTED";
  G4cout << "<<<<<END>>>>>   " << resume << G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmAcceptance::EmAcceptanceGauss(const G4String& title, G4int stat,
                                           G4double avr, G4double avr0,
                                           G4double rms, G4double limit)
{
  G4double x = std::sqrt((G4double)stat);
  G4double dde = avr - avr0;
  G4double de = dde*x/rms;

  G4cout << title << ": " << avr << "  del"<< title << "= " << dde
         << " nrms= " << de << G4endl;

  if (std::fabs(de) > limit) fIsAccepted = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
