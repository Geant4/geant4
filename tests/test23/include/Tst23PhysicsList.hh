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
#ifndef Tst23PhysicsList_h
#define Tst23PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
class Tst23PhysicsListMessenger;
class Tst23PhysicsList: public G4VModularPhysicsList
{
public:
  Tst23PhysicsList();
  virtual ~Tst23PhysicsList();
  
public:
  // SetCuts() 
  virtual void SetCuts();

  void SetCutInRangeForRegion( G4double cut, G4int region);  

private:
  enum { N_Regions      = 2}; 
  void SetCutValue(G4double* cuts, const G4String& name); 
  void SetCutValueForOthers(G4double* cuts);
  void SetParticleCuts( G4double* cuts, G4ParticleDefinition* particle);

private:
  G4double     cutInRange[N_Regions];
  Tst23PhysicsListMessenger* myMessenger;
};


#endif



