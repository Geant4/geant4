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
// Class Description:
// The list of particles and processes are defined in this class.
// Class Description - end
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test17PhysicsList_h
#define Test17PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include <vector>

class Test17DetectorConstruction;
class Test17PhysicsListMessenger;
class Test17StepCut;
class G4hLowEnergyIonisation;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17PhysicsList: public G4VUserPhysicsList
{
public: // Without description

  Test17PhysicsList( Test17DetectorConstruction*);
  virtual ~Test17PhysicsList();

protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();

  virtual void SetCuts();

private:
  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();

public: // Without description

  void SetGammaCut(G4double);
  void SetElectronCut(G4double);
  void SetMaxStep(G4double);
  void SetCutForSecondaryPhotons(G4double);
  void SetCutForAugerElectrons(G4double);

private:

  G4double cutForGamma;
  G4double cutForElectron;

  G4double MaxChargedStep;

  std::vector<G4hLowEnergyIonisation*> hionVector;

  Test17DetectorConstruction* pDet;
  Test17PhysicsListMessenger* physicsListMessenger;
  Test17StepCut* theStepCut;
};

#endif



