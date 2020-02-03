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
/// \file electromagnetic/TestEm10/include/PhysicsList.hh
/// \brief Definition of the PhysicsList class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   PhysicsList
//
// Description: EM physics with a possibility to add transition radiation
//
// Author:      V.Ivanchenko 16.12.2015
//
//----------------------------------------------------------------------------
//

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class PhysicsListMessenger;
class TransitionRadiationPhysics;
class StepMax;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
public:

  PhysicsList(DetectorConstruction* ptr);
  virtual ~PhysicsList();

  virtual void ConstructParticle();

  virtual void ConstructProcess();
    
  void AddPhysicsList(const G4String& name);

  void SetXTRModel(const G4String& name);

private:

  void AddStepMax(); 

  G4VPhysicsConstructor*  fEmPhysicsList;
  G4VPhysicsConstructor*  fDecayPhysicsList;

  TransitionRadiationPhysics* fXTRPhysicsList;

  G4String fEmName;
  G4int    fVerbose;

  PhysicsListMessenger* fMessenger;
  static G4ThreadLocal StepMax* fStepMaxProcess;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

