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
// $Id: PrimaryGeneratorAction.hh,v 1.6 2010-10-26 12:09:14 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  ~PrimaryGeneratorAction();

public:

  //It reads the parameters of the primary particles.
  //Generates the primary event via the ParticleGun method.
  void GeneratePrimaries(G4Event* anEvent);

  //Get/Set methods
  void SetBeamSigmaE(G4double val);
  void SetBeamEnergy(G4double val);

  inline void SetBeamX(G4double val) {x0 = val;};
  inline void SetBeamY(G4double val) {y0 = val;};
  inline void SetBeamZ(G4double val) {z0 = val;};
  inline void SetBeamSigmaX(G4double val) {sigmaX = val;};
  inline void SetBeamSigmaY(G4double val) {sigmaY = val;};
  inline void SetBeamSigmaZ(G4double val) {sigmaY = val;};
  inline void SetBeamMinCosTheta(G4double val) {minCosTheta = val;};
  inline void SetSigmaTheta(G4double val) {sigmaTheta = val;};
  inline void SetVerbose(G4int val) {verbose = val;};
  inline void SetRandom(const G4String& type) {m_gauss = type;};

private:

  void InitializeMe();

  G4ParticleGun* particleGun;
  PrimaryGeneratorMessenger* theMessenger;

  DetectorConstruction* fDetector;

  G4int counter;
  G4int verbose;
  G4double x0, y0, z0;
  G4double sigmaX, sigmaY, sigmaZ;
  G4double rMax2;
  G4double sigmaE;
  G4double sigmaTheta;
  G4double energy;
  G4double minCosTheta;
  G4ThreeVector position;
  G4ThreeVector direction;
  G4String m_gauss;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif






