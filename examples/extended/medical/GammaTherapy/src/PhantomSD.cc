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
// $Id: PhantomSD.cc 103469 2017-04-11 07:29:36Z gcosmo $
//
/// \file medical/GammaTherapy/src/PhantomSD.cc
/// \brief Implementation of the PhantomSD class
//
// -------------------------------------------------------------
//
//
//      ---------- PhantomSD -------------
//              
//  Modified:
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "PhantomSD.hh"
#include "G4RunManager.hh"
#include "Run.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhantomSD::PhantomSD(const G4String& name)
  : G4VSensitiveDetector(name), fShiftZ(0.0),fCounter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhantomSD::~PhantomSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhantomSD::Initialize(G4HCofThisEvent*)
{
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  ++fCounter;
  if(run->GetVerbose()) {
    G4cout << "PhantomSD: Begin Of Event # " << fCounter << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool PhantomSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4double edep = aStep->GetTotalEnergyDeposit();

  // only if there is energy deposition
  if(0.0 < edep) {

    G4ThreeVector p1 = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector p2 = aStep->GetPostStepPoint()->GetPosition();
    G4double x1 = p1.x();
    G4double y1 = p1.y();
    G4double z1 = p1.z() - fShiftZ;
    G4double r1 = std::sqrt(x1*x1 + y1*y1);
    G4double x2 = p2.x();
    G4double y2 = p2.y();
    G4double z2 = p2.z() - fShiftZ;
    G4double r2 = std::sqrt(x2*x2 + y2*y2);
    G4double x0 = 0.5*(x1 + x2);
    G4double y0 = 0.5*(y1 + y2);
    G4double z0 = 0.5*(z1 + z2);
    G4double r0 = std::sqrt(x0*x0 + y0*y0);

    run->AddPhantomStep(edep,r1,z1,r2,z2,r0,z0);

    if(run->GetVerbose()) {
      G4cout << "PhantomSD: energy = " << edep/MeV
             << " MeV is deposited at the step at r1,z1= " << r1 << " " << z1
             << "; r2,z2= " << r2 <<  " " << z2 << G4endl;
    }
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhantomSD::EndOfEvent(G4HCofThisEvent*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhantomSD::clear()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void PhantomSD::PrintAll()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

