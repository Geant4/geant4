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
// $Id: Tst50PrimaryGeneratorAction.hh,v 1.8 2003-05-17 18:11:52 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
 
#ifndef Tst50PrimaryGeneratorAction_h
#define Tst50PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;
class Tst50PrimaryGeneratorMessenger; 

class Tst50PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  Tst50PrimaryGeneratorAction();    
  ~Tst50PrimaryGeneratorAction();
  void GeneratePrimaries(G4Event*);

public:
  G4double GetInitialEnergy();
  G4String GetParticle();
  void SetRandomDirection (G4String val) {randomDirection = val;}
  
private:
  G4String randomDirection;
  G4ParticleGun* particleGun;
  Tst50PrimaryGeneratorMessenger* gunMessenger; 
};
#endif










