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
// $Id: G4ionGasIonisation.cc,v 1.15 2009/11/10 11:50:30 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ionGasIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.07.2007
//
// Modifications:
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ionGasIonisation.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4ionGasIonisation::G4ionGasIonisation(const G4String& name)
  : G4ionIonisation(name)
{
  G4cout << G4endl;
  G4cout << "!!! G4ionGasIonisation class is obsolete and may be removed for the next major Geant4 release !!!" << G4endl;
  G4cout << "!!! Please use G4ionIonisation for ions !!!" << G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionGasIonisation::~G4ionGasIonisation()
{}

/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionGasIonisation::SampleChargeAfterStep(G4double qeff, G4double xeff)
{
  // qeff - equilibrium charge
  // xeff - effective number of collisions
  // q    - current charge
  G4double q = G4double(currentIonZ);
  if(qeff > q) {
    if(G4UniformRand() < qeff - q) currentIonZ++;
  } else {
    if(G4UniformRand() < q - qeff) currentIonZ--;
  }

  q = eplus*currentIonZ;
  if(verboseLevel > 1) G4cout << "G4ionGasIonisation: Q1= " << currentIonZ
			      << " Qeff= " << qeff/eplus << "  Neff= " << xeff
			      << G4endl;
  return q;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
