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
//
// $Id: RemSimPrimaryGeneratorAction.hh,v 1.11 2005-09-08 06:56:18 guatelli Exp $// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli, guatelli@ge.infn.it
//
#ifndef RemSimPrimaryGeneratorAction_h
#define RemSimPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4DataVector.hh"

class G4ParticleGun;
class G4Event;
class RemSimPrimaryGeneratorMessenger;
class RemSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  RemSimPrimaryGeneratorAction();
  ~RemSimPrimaryGeneratorAction();

public:
  G4double GetInitialEnergy();
  void GeneratePrimaries(G4Event* anEvent);
  void SelectPrimaries(G4String value); 
  void Read(G4String);

private:
  void ReadData(G4double, G4String);
  G4double GetPrimaryParticleEnergyDistributionSum();

private:  
  G4ParticleGun* particleGun;
  G4String value;
  RemSimPrimaryGeneratorMessenger* messenger;
  G4DataVector* energies;
  G4DataVector* data;
};
#endif


