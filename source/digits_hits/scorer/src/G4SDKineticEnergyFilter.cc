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
// G4VSensitiveDetector
#include "G4SDKineticEnergyFilter.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"
////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector.
//
//
// Created: 2005-11-14  Tsukasa ASO.
//
///////////////////////////////////////////////////////////////////////////////

G4SDKineticEnergyFilter::G4SDKineticEnergyFilter(G4String name, G4double elow,
                                                 G4double ehigh)
  : G4VSDFilter(name)
  , fLowEnergy(elow)
  , fHighEnergy(ehigh)
{}

G4bool G4SDKineticEnergyFilter::Accept(const G4Step* aStep) const
{
  G4double kinetic = aStep->GetPreStepPoint()->GetKineticEnergy();
  if(kinetic < fLowEnergy)
    return false;
  if(kinetic >= fHighEnergy)
    return false;
  return true;
}

void G4SDKineticEnergyFilter::SetKineticEnergy(G4double elow, G4double ehigh)
{
  fLowEnergy  = elow;
  fHighEnergy = ehigh;
}

void G4SDKineticEnergyFilter::show()
{
  G4cout << " G4SDKineticEnergyFilter:: " << GetName() << " LowE  "
         << G4BestUnit(fLowEnergy, "Energy") << " HighE "
         << G4BestUnit(fHighEnergy, "Energy") << G4endl;
}
