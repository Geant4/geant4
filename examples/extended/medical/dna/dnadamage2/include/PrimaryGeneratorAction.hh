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
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: J. Naoki D. Kondo (UCSF, US)
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef DNADAMAGE2_PrimaryGeneratorAction_h
#define DNADAMAGE2_PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include <G4UImessenger.hh>
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction,
                               public G4UImessenger
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;
  void GeneratePrimaries(G4Event*) override;
  const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
  void SetNewValue(G4UIcommand * command, G4String newValue) override;

private:
  void LoadSpectrum();
  G4double SampleSpectrum();
  G4double SampleProbability(G4double);
  G4ThreeVector SamplePosition();

private:
  G4ParticleGun* fParticleGun = nullptr;
  G4UIcmdWithAString* fpSourceFileUI = nullptr;
  G4UIcmdWithAnInteger* fpVertexUI  = nullptr;
  G4String fSourceFile = "";
  G4int fVertex = 1;
  std::vector<G4double> fEnergyWeights;
  std::vector<G4double> fEnergies;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
