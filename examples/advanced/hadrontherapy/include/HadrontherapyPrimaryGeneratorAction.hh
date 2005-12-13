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
// $Id: HadrontherapyPrimaryGeneratorAction.hh; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

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
 

  void SetsigmaEnergy(G4double);
  void SetmeanKineticEnergy(G4double);
  void SetDefaultPrimaryParticle();
  void GeneratePrimaries(G4Event*);
  void SetXposition(G4double);
  void SetYposition(G4double);
  void SetZposition(G4double);
  void SetsigmaY(G4double);
  void SetsigmaZ(G4double);
  void SetsigmaMomentumY(G4double);
  void SetsigmaMomentumZ(G4double);
  
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
  G4ParticleGun*                particleGun;
  HadrontherapyPrimaryGeneratorMessenger* gunMessenger; 
  G4double sigmaX;
};

#endif


