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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
//

#ifndef GammaRayTelPhysicsList_h
#define GammarayTelPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4EmConfigurator.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class GammaRayTelPhysicsListMessenger;

class GammaRayTelPhysicsList: public G4VModularPhysicsList {
public:
	explicit GammaRayTelPhysicsList();

	~GammaRayTelPhysicsList() override;

	void AddPackage(const G4String &name);

	void AddPhysicsList(const G4String &name);

	void ConstructParticle() override;

	void ConstructProcess() override;

	void SetCutForElectron(G4double cut);

	void SetCutForGamma(G4double cut);

	void SetCutForPositron(G4double cut);

private:
	G4EmConfigurator emConfigurator;

	G4bool helIsRegisted;

	G4bool bicIsRegisted;

	G4bool biciIsRegisted;

	G4bool locIonIonInelasticIsRegistered;

	G4bool radioactiveDecayIsRegisted;

	G4String emName;

	G4VPhysicsConstructor *emPhysicsList;

	G4VPhysicsConstructor *decPhysicsList;

	std::vector<G4VPhysicsConstructor*> hadronPhys;

	GammaRayTelPhysicsListMessenger *pMessenger;
};
#endif
