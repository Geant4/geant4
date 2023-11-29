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
/// \file B5/include/PrimaryGeneratorAction.hh
/// \brief Definition of the B5::PrimaryGeneratorAction class

#ifndef B5PrimaryGeneratorAction_h
#define B5PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <CLHEP/Units/SystemOfUnits.h>

class G4ParticleGun;
class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

namespace B5
{

/// Primary generator
///
/// A single particle is generated.
/// User can select
/// - the initial momentum and angle
/// - the momentum and angle spreads
/// - random selection of a particle type from proton, kaon+, pi+, muon+, e+


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

    void SetMomentum(G4double val) { fMomentum = val; }
    G4double GetMomentum() const { return fMomentum; }

    void SetSigmaMomentum(G4double val) { fSigmaMomentum = val; }
    G4double GetSigmaMomentum() const { return fSigmaMomentum; }

    void SetSigmaAngle(G4double val) { fSigmaAngle = val; }
    G4double GetSigmaAngle() const { return fSigmaAngle; }

    void SetRandomize(G4bool val) { fRandomizePrimary = val; }
    G4bool GetRandomize() const { return fRandomizePrimary; }

  private:
    void DefineCommands();

    G4ParticleGun* fParticleGun = nullptr;
    G4GenericMessenger* fMessenger = nullptr;
    G4ParticleDefinition* fPositron = nullptr;
    G4ParticleDefinition* fMuon = nullptr;
    G4ParticleDefinition* fPion = nullptr;
    G4ParticleDefinition* fKaon = nullptr;
    G4ParticleDefinition* fProton = nullptr;
    G4double fMomentum = 1000. * CLHEP::MeV;
    G4double fSigmaMomentum = 50. * CLHEP::MeV;
    G4double fSigmaAngle = 2. * CLHEP::deg;
    G4bool fRandomizePrimary = true;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
