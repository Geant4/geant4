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
//
// $Id: PhotInPrimaryGeneratorAction.hh,v 1.5 2006/06/29 16:24:53 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//

#ifndef PhotInPrimaryGeneratorAction_h
#define PhotInPrimaryGeneratorAction_h 1

#include "PhotInDetectorConstruction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

class PhotInPrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction
{
public:
  PhotInPrimaryGeneratorAction();    
  virtual ~PhotInPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

  //void SetSerial(G4bool ser) {serial  = ser;} // Different positions for different setups
  void SetDetector(PhotInDetectorConstruction* det);
  void SetProjectileName(G4String partName);
  void SetProjectileEnergy(G4double partEnergy);
  void SetSection(G4int sec)                  // Define position for the particular section
  {
    if(sec<0||sec>PhotInNumSections)
    {
      G4cout<<"PhotInPrimaryGeneratorAction::SetSection: section="<<sec<<"? set 1"<<G4endl;
      section = 1;
    }
				else section = sec;
  }
  //G4bool GetSerial() {return serial;} // Get the setups
  G4int GetSection() {return section;} // Get the starting section
  G4String GetProjectileName() {return part;} // Get the projectile name
  G4double GetProjectileEnergy() {return energy;} // Get the projectile energy

private: //--- BODY ---
  G4ParticleGun*                particleGun;
  G4int                         section;
  PhotInDetectorConstruction*   detector;
  G4String                      part;
  G4double                      energy;
  G4String                      oldPart;
  G4double                      oldEnergy;
};

#endif


