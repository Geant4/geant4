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
// This class was made using G4MoleculeGun as a base
// G4MoleculeGun Author: Mathieu Karamitros
//
// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508
//
/// \file MoleculeInserter.hh
/// \brief Definition of the DNA chemical species inserter for IRT

#ifndef DNADAMAGE2_moleculeinserter_h
#define DNADAMAGE2_moleculeinserter_h 1

#include "G4ITGun.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <map>
#include <G4memory.hh>

class G4Track;
class MoleculeInserter;
class G4Molecule;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MoleculeShoot {
public:
  MoleculeShoot();
  ~MoleculeShoot() = default;
  void Shoot(MoleculeInserter*);

public:
  G4int    fNumber = 0;
  G4String fMoleculeName = "None";
  G4double fTime = 1;
  G4ThreeVector fPosition = G4ThreeVector();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MoleculeInserter : public G4ITGun {
public:
  MoleculeInserter();
  MoleculeInserter(G4bool);
  ~MoleculeInserter() = default;
  void DefineTracks() override;
  void AddMolecule(const G4String& moleculeName, 
                   const G4ThreeVector& position, double time = 0);
  void PushToChemistry(G4Molecule*, G4double, G4ThreeVector);
  void CreateMolecule(G4Molecule*, G4double, G4ThreeVector);
  void CreateMolecule(G4String, G4double, G4ThreeVector);
  void Clean();

  std::vector<G4Track*> const GetInsertedTracks() {return fInsertedTracks;}

private:
  G4bool fSaveTrackID = true;
  std::vector<G4Track*> fInsertedTracks;

protected:
  std::vector<MoleculeShoot> fShoots;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
