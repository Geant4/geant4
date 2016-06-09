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
//
// $Id: EmAcceptance.cc,v 1.3 2004/12/02 19:06:05 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EmAcceptance.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmAcceptance::EmAcceptance()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EmAcceptance::~EmAcceptance()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmAcceptance::BeginOfAcceptance(const G4String& title, G4int stat)
{
  G4cout << G4endl;
  G4cout << "<<<<<ACCEPTANCE>>>>> " << stat << " events for " << title << G4endl;
  isAccepted = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EmAcceptance::EndOfAcceptance()
{
  G4String resume = "IS ACCEPTED";
  if(!isAccepted) resume = "IS NOT ACCEPTED";
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

  G4cout << title << ": " << avr << "  del"<< title << "= " << dde << " nrms= " << de << G4endl;
  if(std::fabs(de) > limit) isAccepted = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
