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
// --------------------------------------------------------------
// Authors: E. Delage
// november 2013
// --------------------------------------------------------------
//
//
/// \file DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

/// Messenger class that defines commands for DetectorConstruction.
/// - /PDB4DNA/det/loadPDB stringfilename
/// - /PDB4DNA/det/drawAtoms
/// - /PDB4DNA/det/drawNucleotides
/// - /PDB4DNA/det/drawResidues
/// - /PDB4DNA/det/buildBoundingV
/// - /PDB4DNA/det/drawAtomsWithBounding
/// - /PDB4DNA/det/drawNucleotidesWithBounding
/// - /PDB4DNA/det/drawResiduesWithBounding

class DetectorMessenger: public G4UImessenger
{
public:
  DetectorMessenger(DetectorConstruction*);
  virtual ~DetectorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

private:
  DetectorConstruction*  fpDetectorConstruction;

  G4UIdirectory*                fpDirectory;
  G4UIdirectory*                fpDetDirectory;
  G4UIcmdWithAString*           fpLoadPdbCmd;
  G4UIcmdWithoutParameter*      fpBuildBoundingV;

  G4UIcmdWithoutParameter*      fpDrawAtoms;
  G4UIcmdWithoutParameter*      fpDrawNucleotides;
  G4UIcmdWithoutParameter*      fpDrawResidues;
  G4UIcmdWithoutParameter*      fpDrawAtomsWithBounding;
  G4UIcmdWithoutParameter*      fpDrawNucleotidesWithBounding;
  G4UIcmdWithoutParameter*      fpDrawResiduesWithBounding;
};

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
