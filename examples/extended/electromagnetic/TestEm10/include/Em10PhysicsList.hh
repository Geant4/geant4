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
/// \file electromagnetic/TestEm10/include/Em10PhysicsList.hh
/// \brief Definition of the Em10PhysicsList class
//
//
// $Id: Em10PhysicsList.hh 66241 2012-12-13 18:34:42Z gunter $
//

#ifndef Em10PhysicsList_h
#define Em10PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4ForwardXrayTR ;
class Em10StepCut;
class Em10DetectorConstruction;
class Em10PhysicsListMessenger;
class G4ProductionCuts;


class Em10PhysicsList: public G4VModularPhysicsList  // G4VUserPhysicsList
{
public:
  Em10PhysicsList( Em10DetectorConstruction*);

  ~Em10PhysicsList();

  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
 
  void SetCuts();

private:
    // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBarions();

  // these methods Construct physics processes and register them

  void AddParameterisation();
  void ConstructGeneral();
  void ConstructEM();

public:

  void SetGammaCut(G4double);
  void SetElectronCut(G4double);

  void SetRegGammaCut(G4double    cut){fGammaCut    = cut;};
  void SetRegElectronCut(G4double cut){fElectronCut = cut;};
  void SetRegPositronCut(G4double cut){fPositronCut = cut;};

  void SetRadiatorCuts();
  void SetDetectorCuts();

  void SetMaxStep(G4double);
  void SetXTRModel(G4String m) {fXTRModel = m; G4cout<<fXTRModel<<G4endl;}; 

private:

  G4double MaxChargedStep;

  G4ForwardXrayTR*       fForwardXrayTR ;

  Em10StepCut* theeminusStepCut ;
  Em10StepCut* theeplusStepCut ;

  G4double cutForGamma;
  G4double cutForElectron, cutForPositron;

  Em10DetectorConstruction* pDet;

  Em10PhysicsListMessenger* physicsListMessenger;

  G4ProductionCuts* fRadiatorCuts;
  G4ProductionCuts* fDetectorCuts;
  G4double fElectronCut, fGammaCut, fPositronCut;
  G4String fXTRModel;
};

#endif



