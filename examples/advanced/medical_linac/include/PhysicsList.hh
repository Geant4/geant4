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
// Code developed by
// Silvia Pozzi (1), silvia.pozzi@iss.it
// Barbara Caccia (1), barbara.caccia@iss.it
// Carlo Mancini Terracciano (2), carlo.mancini.terracciano@roma1.infn.it
// (1) Istituto Superiore di Sanita' and INFN Roma, Italy
// (2) Univ. La Sapienza and INFN Roma, Italy

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include <G4VModularPhysicsList.hh>
#include "globals.hh"

class G4VPhysicsConstructor;
class StepMax;
class PhysicsListMessenger;

class PhysicsList : public G4VModularPhysicsList
{
public:
	PhysicsList();
	~PhysicsList();

	void SetCuts () override;
	void ConstructParticle() override;
	void AddPhysicsList(const G4String& name);
	void ConstructProcess() override;

	void AddStepMax();


private:

  G4bool helIsRegistered;
  G4bool bicIsRegistered;
  G4bool biciIsRegistered;
  G4bool locIonIonInelasticIsRegistered;

  G4String                             emName;
  G4VPhysicsConstructor*               emPhysicsList;
  G4VPhysicsConstructor*               decPhysicsList;
  std::vector<G4VPhysicsConstructor*>  hadronPhys;

  StepMax* stepMaxProcess;

  PhysicsListMessenger* pMessenger;
};

#endif
