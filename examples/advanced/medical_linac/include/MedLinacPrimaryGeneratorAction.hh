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
// $Id: MedLinacPrimaryGeneratorAction.hh,v 1.3 2006/06/29 16:04:01 gunter Exp $
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


