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
// $Id: PrimaryGeneratorAction.hh 103662 2017-04-20 14:58:33Z gcosmo $
//
/// \file medical/GammaTherapy/include/PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class
//

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

//---------------------------------------------------------------------------
//
// ClassName:   PrimaryGeneratorAction
//
// Description: Generate primary beam
//
// Authors: V.Grichine, V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PrimaryGeneratorMessenger;
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  // The constructor defines a ParticleGun object, which allows
  // shooting a beam of particles through the experimental set-up.
  PrimaryGeneratorAction(DetectorConstruction* pDet);

  //The destructor. It deletes the ParticleGun.
  virtual ~PrimaryGeneratorAction();

  //It reads the parameters of the primary particles.
  //Generates the primary event via the ParticleGun method.
  void GeneratePrimaries(G4Event* anEvent);

  //Get/Set methods
  void SetBeamEnergy(G4double val);

  inline void SetBeamSigmaE(G4double val) { fSigmaE = val; };
  inline void SetBeamX(G4double val) { fX0 = val;};
  inline void SetBeamY(G4double val) { fY0 = val;};
  inline void SetBeamZ(G4double val) { fZ0 = val;};
  inline void SetBeamSigmaX(G4double val) { fSigmaX = val;};
  inline void SetBeamSigmaY(G4double val) { fSigmaY = val;};
  inline void SetBeamSigmaZ(G4double val) { fSigmaY = val;};
  inline void SetBeamMinCosTheta(G4double val) { fMinCosTheta = val;};
  inline void SetSigmaTheta(G4double val) { fSigmaTheta = val;};
  inline void SetVerbose(G4int val) { fVerbose = val;};
  inline void SetRandom(const G4String& type) { fGauss = type;};

  G4bool GetVerbose() const { return fVerbose; }

private:

  void InitializeMe();

  PrimaryGeneratorAction & operator=(const PrimaryGeneratorAction &right);
  PrimaryGeneratorAction(const PrimaryGeneratorAction&);

  G4int fVerbose;
  PrimaryGeneratorMessenger* fMessenger;
  G4ParticleGun* fParticleGun;
  DetectorConstruction* fDetector;

  G4int fCounter;
  G4double fX0, fY0, fZ0;
  G4double fSigmaX, fSigmaY, fSigmaZ;
  G4double fRMax2;
  G4double fSigmaE;
  G4double fSigmaTheta;
  G4double fEnergy;
  G4double fMinCosTheta;
  G4ThreeVector fPosition;
  G4ThreeVector fDirection;
  G4String fGauss;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif






