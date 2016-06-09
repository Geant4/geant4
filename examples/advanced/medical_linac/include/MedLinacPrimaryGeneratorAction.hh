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
// $Id: MedLinacPrimaryGeneratorAction.hh,v 1.2 2004/04/02 17:48:41 mpiergen Exp $
//
//
// Code developed by: M. Piergentili

#ifndef MedLinacPrimaryGeneratorAction_h
#define MedLinacPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <vector>

class G4ParticleGun;
class G4Run;
class G4Event;
class MedLinacDetectorConstruction;
class MedLinacPrimaryGeneratorMessenger;
class MedLinacAnalysisManager;

class MedLinacPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    MedLinacPrimaryGeneratorAction();
    ~MedLinacPrimaryGeneratorAction();

  public:
  void GeneratePrimaries(G4Event* anEvent);
  void SetEnergy(G4double val){ pEnergy = val;}
  void SetSourceType(G4double val) { sigma = val;}
  private:
  G4ParticleGun*                particleGun;  //pointer a to G4  class
  MedLinacPrimaryGeneratorMessenger* gunMessenger;//messenger of this class
  G4double                      pEnergy;
  G4double                      sigma;
};

#endif


