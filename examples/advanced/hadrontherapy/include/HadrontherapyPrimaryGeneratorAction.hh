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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#ifndef HadrontherapyPrimaryGeneratorAction_h
#define HadrontherapyPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class HadrontherapyPrimaryGeneratorMessenger;
class HadrontherapyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  HadrontherapyPrimaryGeneratorAction();    
  ~HadrontherapyPrimaryGeneratorAction();
  
public:
  // Methods to change the parameters of primary particle generation 
  // interactively
  void SetsigmaEnergy(G4double);
  void SetmeanKineticEnergy(G4double);
  void GeneratePrimaries(G4Event*);
  void SetXposition(G4double);
  void SetYposition(G4double);
  void SetZposition(G4double);
  void SetsigmaY(G4double);
  void SetsigmaZ(G4double);
  void SetsigmaMomentumY(G4double);
  void SetsigmaMomentumZ(G4double);
  G4double GetmeanKineticEnergy(void);
  G4ParticleGun *GetParticleGun(void){return particleGun;}
    
private:
  void SetDefaultPrimaryParticle();
  G4double meanKineticEnergy;
  G4double sigmaEnergy;
  G4double X0;
  G4double Y0;
  G4double Z0;
  G4double sigmaY;
  G4double sigmaZ;
  G4double sigmaMomentumY;
  G4double sigmaMomentumZ;

private:
  G4ParticleGun*    		          particleGun;
  HadrontherapyPrimaryGeneratorMessenger* gunMessenger; 
  G4double				  sigmaX;
};

#endif


