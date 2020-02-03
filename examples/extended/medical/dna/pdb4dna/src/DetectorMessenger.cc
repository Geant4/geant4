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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
: G4UImessenger(),
  fpDetectorConstruction(Det)
{ 
  fpDirectory = new G4UIdirectory("/PDB4DNA/");
  fpDirectory->SetGuidance("UI commands of this example");

  fpDetDirectory = new G4UIdirectory("/PDB4DNA/det/");
  fpDetDirectory->SetGuidance("Detector construction control");

  fpLoadPdbCmd = new G4UIcmdWithAString("/PDB4DNA/det/loadPDB",this);
  fpLoadPdbCmd->SetGuidance("Load PDB file");
  fpLoadPdbCmd->SetParameterName("filename",false);
  fpLoadPdbCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fpDrawAtoms = new G4UIcmdWithoutParameter("/PDB4DNA/det/drawAtoms",this);
  fpDrawAtoms->SetGuidance("Draw atoms");
  fpDrawAtoms->AvailableForStates(G4State_PreInit, G4State_Idle);

  fpDrawNucleotides = new G4UIcmdWithoutParameter(
      "/PDB4DNA/det/drawNucleotides",
      this);
  fpDrawNucleotides->SetGuidance("Draw nucleotides with bounding sphere");
  fpDrawNucleotides->AvailableForStates(G4State_PreInit, G4State_Idle);

  fpDrawResidues = new G4UIcmdWithoutParameter(
      "/PDB4DNA/det/drawResidues",
      this);
  fpDrawResidues->SetGuidance("Draw residues inside nucleotides with sphere "
      "linked by cylinders");
  fpDrawResidues->AvailableForStates(G4State_PreInit, G4State_Idle);

  fpBuildBoundingV = new G4UIcmdWithoutParameter(
      "/PDB4DNA/det/buildBoundingV",
      this);
  fpBuildBoundingV->SetGuidance("Build molecule bounding volume");
  fpBuildBoundingV->AvailableForStates(G4State_PreInit, G4State_Idle);

  fpDrawAtomsWithBounding = new G4UIcmdWithoutParameter(
      "/PDB4DNA/det/drawAtomsWithBounding",
      this);
  fpDrawAtomsWithBounding->SetGuidance("Draw atoms with bounding volume");
  fpDrawAtomsWithBounding->AvailableForStates(G4State_PreInit, G4State_Idle);

  fpDrawNucleotidesWithBounding = new G4UIcmdWithoutParameter(
      "/PDB4DNA/det/drawNucleotidesWithBounding",
      this);
  fpDrawNucleotidesWithBounding->SetGuidance(
      "Draw nucleotides with bounding sphere and bounding volume");
  fpDrawNucleotidesWithBounding->AvailableForStates(G4State_PreInit,
                                                    G4State_Idle);

  fpDrawResiduesWithBounding = new G4UIcmdWithoutParameter(
      "/PDB4DNA/det/drawResiduesWithBounding",
      this);
  fpDrawResiduesWithBounding->SetGuidance("Draw residues inside nucleotides"
      " with sphere linked by cylinders and with bounding volume");
  fpDrawResiduesWithBounding->AvailableForStates(
      G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fpLoadPdbCmd;
  delete fpDrawAtoms;
  delete fpDrawNucleotides;
  delete fpDrawResidues;
  delete fpBuildBoundingV;
  delete fpDrawAtomsWithBounding;
  delete fpDrawNucleotidesWithBounding;
  delete fpDrawResiduesWithBounding;

  delete fpDetDirectory;
  delete fpDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fpLoadPdbCmd )
  {
    fpDetectorConstruction->LoadPDBfile(newValue);
  }
  if( command == fpDrawAtoms )
  {
    fpDetectorConstruction->DrawAtoms_();
  }
  if( command == fpDrawNucleotides )
  {
    fpDetectorConstruction->DrawNucleotides_();
  }
  if( command == fpDrawResidues )
  {
    fpDetectorConstruction->DrawResidues_();
  }
  if( command == fpBuildBoundingV )
  {
    fpDetectorConstruction->BuildBoundingVolume();
  }
  if( command == fpDrawAtomsWithBounding )
  {
    fpDetectorConstruction->DrawAtomsWithBounding_();
  }
  if( command == fpDrawNucleotidesWithBounding )
  {
    fpDetectorConstruction->DrawNucleotidesWithBounding_();
  }
  if( command == fpDrawResiduesWithBounding )
  {
    fpDetectorConstruction->DrawResiduesWithBounding_();
  }
}
