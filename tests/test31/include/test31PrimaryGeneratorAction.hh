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
#ifndef test31PrimaryGeneratorAction_h
#define test31PrimaryGeneratorAction_h 1

//---------------------------------------------------------------------------
//
// ClassName:   test31PrimaryGeneratorAction
//  
// Description: Generate primary beam 
//
// Authors:    0.6.04.01 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class test31PrimaryGeneratorMessenger;

class test31PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

  // The constructor defines a ParticleGun object, which allows 
  // shooting a beam of particles through the experimental set-up.
    test31PrimaryGeneratorAction(test31DetectorConstruction*);

  //The destructor. It deletes the ParticleGun.
    ~test31PrimaryGeneratorAction();

  public:

  //It reads the parameters of the primary particles.
  //Generates the primary event via the ParticleGun method.
    void GeneratePrimaries(G4Event* anEvent);

  //Get/Set methods
    void SetBeamX(G4double val) {x0 = val;};
    void SetBeamY(G4double val) {y0 = val;};
    void SetBeamZ(G4double val) {z0 = val;};
    void SetBeamSigmaX(G4double val) {sigmaX = val;};
    void SetBeamSigmaY(G4double val) {sigmaY = val;};
    void SetBeamSigmaZ(G4double val) {sigmaY = val;};
    void SetBeamSigmaE(G4double val); 
    void SetBeamEnergy(G4double val); 
    void SetBeamMinCosTheta(G4double val) {minCosTheta = val;};
    void SetVerbose(G4int val) {verbose = val;};
    G4ThreeVector GetBeamPosition() const {return position;}; 
    G4ThreeVector GetBeamDirection() const {return direction;};
    G4ThreeVector GetBeamEnergy() const {return energy;};
    G4int GetVerbose() const {return verbose;};
    void SetBeamBeta(G4double val);
    void SetSigmaBeta(G4double val);
    void SetRandom(const G4String& type) {m_gauss = type;};

  private:

    void InitializeMe();

    test31DetectorConstruction* theDet;
    G4ParticleGun* particleGun;
    test31PrimaryGeneratorMessenger* theMessenger;

    G4int counter;
    G4int verbose;
    G4double x0, y0, z0;
    G4double sigmaX, sigmaY, sigmaZ;
    G4double sigmaE;
    G4double energy;
    G4double minE, maxE, minBeta, maxBeta;
    G4double minCosTheta;
    G4ThreeVector position;
    G4ThreeVector direction;
    G4String m_gauss;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif






