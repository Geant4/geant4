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
// MoleculeInserter.cc
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file MoleculeInserter.cc
/// \brief Implementation of the DNA chemical species inserter for IRT

#include "MoleculeInserter.hh"
#include "G4MoleculeTable.hh"
#include "G4Molecule.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include <cassert>
#include "Randomize.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Track.hh"
#include "G4Molecule.hh"
#include "G4ITTrackHolder.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

MoleculeInserter::MoleculeInserter():fSaveTrackID(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

MoleculeInserter::MoleculeInserter(G4bool save):fSaveTrackID(save)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void MoleculeShoot::Shoot(MoleculeInserter* gun) 
{
  for(int i = 0; i < fNumber; ++i) {
    G4MolecularConfiguration* conf = 
      G4MoleculeTable::Instance()->GetConfiguration(fMoleculeName);

    if (conf == 0) {
      G4String msg = "Chemistry Error: Molecule " + fMoleculeName + " don't exists.";
      G4Exception("MoleculeInserter::Shoot()", "Invalid_Value", FatalException, msg);
    }

    G4Molecule* molecule = new G4Molecule(conf);
    gun->PushToChemistry(molecule, fTime, fPosition);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void MoleculeInserter::CreateMolecule(G4Molecule* molecule, 
                                      G4double time, G4ThreeVector pos)
{
  G4Track* MolTrack = molecule->BuildTrack(time, pos);
  PushTrack(MolTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void MoleculeInserter::CreateMolecule(G4String molecule,
                                      G4double time, G4ThreeVector pos)
{
  G4MolecularConfiguration* conf = 
          G4MoleculeTable::Instance()->GetConfiguration(molecule);

  if (conf == 0) {
    G4String msg = "Chemistry Error: Molecule " + molecule + " don't exists.";
    G4Exception("MoleculeInserter::CreateMolecule()", 
                "Invalid_Value", FatalException, msg);
  }

  G4Molecule* gmolecule = new G4Molecule(conf);
  G4Track* MolTrack = gmolecule->BuildTrack(time, pos);
  PushTrack(MolTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void MoleculeInserter::PushToChemistry(G4Molecule* mol,
                                       G4double time, G4ThreeVector position)
{
  G4Track* MolTrack = mol->BuildTrack(time, position);
  PushTrack(MolTrack);

  if (fSaveTrackID) {
    fInsertedTracks.push_back(MolTrack);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void MoleculeInserter::DefineTracks() 
{
  if (fInsertedTracks.size() != 0) {fInsertedTracks.clear();}
  
  for (size_t i = 0; i < fShoots.size(); i++) {
    fShoots[i].Shoot(this);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void MoleculeInserter::AddMolecule(const G4String& name, 
             const G4ThreeVector& position, double time)
{
  MoleculeShoot shoot;
  shoot.fMoleculeName = name;
  shoot.fPosition     = position;
  shoot.fTime         = time;
  fShoots.push_back(shoot);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

MoleculeShoot::MoleculeShoot() 
{
  fMoleculeName = "";
  fTime = 0;
  fNumber = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void MoleculeInserter::Clean() 
{
  if (fShoots.size() != 0)
    fShoots.clear();

  if (fInsertedTracks.size() != 0)
    fInsertedTracks.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
