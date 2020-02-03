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
//  Author: F. Poignant, floriane.poignant@gmail.com
//



#ifndef STCyclotronPhysicsList_h
#define STCyclotronPhysicsList_h 1

#include "STCyclotronDetectorConstruction.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include <vector>

class G4PhysicsConstructor;
class STCyclotronDetectorConstruction;

class STCyclotronPhysicsList: public G4VModularPhysicsList
{
public:
 
  STCyclotronPhysicsList(STCyclotronDetectorConstruction* det);
  virtual ~STCyclotronPhysicsList();
  void ConstructParticle();
  void ConstructProcess();
  
  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
  void SetCutFoil(G4double cutProton, G4double cutElectron, G4double cutPositron, G4double cutGamma, G4double cutNeutron);
  void SetCutTarget(G4double cutProton, G4double cutElectron, G4double cutPositron, G4double cutGamma, G4double cutNeutron);
  
 


  
private:
  G4String                   fEmName;
  G4VPhysicsConstructor*     fEmPhysicsList;
  G4VPhysicsConstructor*     fDecPhysicsList;
  G4VPhysicsConstructor*     fHadPhysicsList;
  G4VPhysicsConstructor*     fRaddecayList;
 
  G4double fThickness_target;
  G4double fThickness_foil;

  G4double fCutForGamma;
  G4double fCutForElectron;
  G4double fCutForPositron;
  
  G4double fCutTargetProton;
  G4double fCutTargetElectron;
  G4double fCutTargetPositron;
  G4double fCutTargetGamma;
  G4double fCutTargetNeutron;

  G4double fCutFoilProton;
  G4double fCutFoilElectron;
  G4double fCutFoilPositron;
  G4double fCutFoilGamma;
  G4double fCutFoilNeutron;
  
  STCyclotronDetectorConstruction* fDetector;

};

#endif
