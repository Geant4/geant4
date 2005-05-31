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
// $Id: PhotInPrimaryGeneratorAction.hh,v 1.2 2005-05-31 15:23:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

private: //--- BODY ---
  G4ParticleGun*                particleGun;
  G4int                         section;
  PhotInDetectorConstruction*   detector;

};

#endif


