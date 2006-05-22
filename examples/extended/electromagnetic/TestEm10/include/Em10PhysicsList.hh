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
// $Id: Em10PhysicsList.hh,v 1.10 2006-05-22 19:05:49 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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



