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
#ifndef test31VEMPhysicsList_h
#define test31VEMPhysicsList_h 1

//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
//
// ClassName:   test31VEMPhysicsList
//  
// Description: Virtual class to build Hadron Physics List for Geant4
//
// Authors:    08.04.01 V.Ivanchenko 
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VUserPhysicsList.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class test31VEMPhysicsList:  public G4VUserPhysicsList
{
public:
  test31VEMPhysicsList() {};
  ~test31VEMPhysicsList() {};

public:
  void ConstructEM() {ConstructProcess();};
  void ConstructParticle() {};
  void SetCuts() {};

  inline void SetVerbose(G4int val) {verbose = val;};
  inline void SetNuclearStopping(G4bool val) {nuclStop = val;};
  inline void SetBarkas(G4bool val) {barkas = val;};
  inline void SetEStoppingTable(G4String name) {table = name;};
  
private:

  // hide assignment operator 
  test31VEMPhysicsList & operator=(const test31VEMPhysicsList &right);
  test31VEMPhysicsList(const test31VEMPhysicsList&);

protected:

  G4int verbose;
  G4bool nuclStop;
  G4bool barkas;
  G4String table;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


