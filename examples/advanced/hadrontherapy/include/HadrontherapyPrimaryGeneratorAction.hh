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
// $Id: HadrontherapyPrimaryGeneratorAction.hh,v 1.0
//
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// -----------------------------------------------------------------

#ifndef HadrontherapyPrimaryGeneratorAction_h
#define HadrontherapyPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class HadrontherapyDetectorConstruction;
class HadrontherapyPrimaryGeneratorMessenger;
class HadrontherapyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction

// ---------------------------------------------------------------------------
{
public:
  HadrontherapyPrimaryGeneratorAction(HadrontherapyDetectorConstruction*);    
  ~HadrontherapyPrimaryGeneratorAction();
  
public:
  void SetDefaultKinematic();
  void GeneratePrimaries(G4Event*);

private:
  G4ParticleGun*                particleGun;
  HadrontherapyDetectorConstruction*      HadrontherapyDetector;
  
  G4double X0;
  G4double Y0;
  G4double Z0;

  G4double sigmaX;
  G4double sigmaY;
  G4double sigmaZ;

  HadrontherapyPrimaryGeneratorMessenger* gunMessenger;     
};

#endif


