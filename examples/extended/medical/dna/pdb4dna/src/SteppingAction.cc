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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "G4EventManager.hh"
#include "EventAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
:G4UserSteppingAction(),RunInitObserver(),fpEventAction(0),fpDetector(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::Initialize()
{
  fpEventAction = (EventAction*) G4EventManager::GetEventManager()->
      GetUserEventAction();
  fpDetector = (DetectorConstruction*)G4RunManager::GetRunManager()->
      GetUserDetectorConstruction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* theStep)
{
  if(theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!=
      "Transportation")
  {
    // Get position and edep of current step
    //
    G4double x = theStep->GetPreStepPoint()->GetPosition().x()/nanometer;
    G4double y = theStep->GetPreStepPoint()->GetPosition().y()/nanometer;
    G4double z = theStep->GetPreStepPoint()->GetPosition().z()/nanometer;
    G4double edepStep = theStep->GetTotalEnergyDeposit()/eV;

    G4LogicalVolume* targetVolume =
        G4LogicalVolumeStore::GetInstance()->GetVolume("BoundingLV");
    G4LogicalVolume* theVolume =
        theStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();

    if ((edepStep > 0.) && (theVolume==targetVolume))
    {
      // Add edep to this event
      //
      fpEventAction->AddEdepEvent(edepStep);
      if (fpDetector->GetBarycenterList()==NULL)
      {
        G4cout << "Barycenter list is null!!!" << G4endl;
      }
      else
      {
        CheckAndProcessDNAHit(x,y,z,edepStep);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SteppingAction::CheckAndProcessDNAHit(G4double x,G4double y, G4double z,
    G4double edepStep)
{
  int numStrand=0;
  int numNucl=0;
  int intResidue=-1; // 0 for Phospat, 1 for Sugar, 2 for Base
  unsigned short int hit = (fpDetector->GetPDBlib()).ComputeMatchEdepDNA(
      fpDetector->GetBarycenterList(),
      fpDetector->GetMoleculeList(),
      x*10., y*10., z*10.,// x10 => angstrom<->nm
      numStrand, numNucl, intResidue);

  if (hit==1)
  {
    if ((intResidue==0)||(intResidue==1)) //Edep in Phosphate or Sugar
    {
      fpEventAction->AddEdepToNucleotide(numStrand,numNucl,edepStep);
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}
