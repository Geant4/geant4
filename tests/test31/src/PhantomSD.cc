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
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Positron.hh"
#include "globals.hh"
#include "Histo.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhantomSD::PhantomSD(const G4String& name)
 :G4VSensitiveDetector(name),
  theHisto(Histo::GetPointer()),
  evno(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhantomSD::~PhantomSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhantomSD::Initialize(G4HCofThisEvent*)
{
  evno++;
  if(0 < theHisto->GetVerbose())
    G4cout << "PhantomSD: Begin Of Event # " << evno << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool PhantomSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  theHisto->AddStep();
  if(0.0 == edep) return true;

  G4ThreeVector p1 = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector p2 = aStep->GetPostStepPoint()->GetPosition();
  G4double x1 = p1.x();
  G4double y1 = p1.y();
  G4double z1 = p1.z() - shiftZ;
  G4double r1 = std::sqrt(x1*x1 + y1*y1);
  G4double x2 = p2.x();
  G4double y2 = p2.y();
  G4double z2 = p2.z() - shiftZ;
  G4double r2 = std::sqrt(x2*x2 + y2*y2);
  G4double x0 = 0.5*(x1 + x2);
  G4double y0 = 0.5*(y1 + y2);
  G4double z0 = 0.5*(z1 + z2);
  G4double r0 = std::sqrt(x0*x0 + y0*y0);

  theHisto->AddStep(edep,r1,z1,r2,z2,r0,z0);

  if(1 < theHisto->GetVerbose()) {
      G4cout << "PhantomSD: energy = " << edep/MeV
             << " MeV is deposited at the step at r1,z1= " << r1 << " " << z1
             << "; r2,z2= " << r2 <<  " " << z2 << G4endl;
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

