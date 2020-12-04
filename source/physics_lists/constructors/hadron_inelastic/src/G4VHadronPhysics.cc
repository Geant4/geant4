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
//---------------------------------------------------------------------------
//
// ClassName:  G4VHadronPhysics
//
// Author: 28 June 2009 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4VHadronPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicProcessType.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronElasticProcess.hh"
#include "G4Neutron.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4VHadronModelBuilder.hh"

G4VHadronPhysics::G4VHadronPhysics(const G4String& aName, G4int verb)
  : G4VPhysicsConstructor(aName)
{
  SetVerboseLevel(verb);
  if (verboseLevel>1) {
    G4cout << "### G4VHadronPhysics: <" << aName << "> is created "
	   << G4endl;
  }
}

G4VHadronPhysics::~G4VHadronPhysics() 
{}

void G4VHadronPhysics::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}

G4HadronicInteraction* 
G4VHadronPhysics::BuildModel(G4VHadronModelBuilder* mBuilder,
			     G4double emin, 
			     G4double emax)
{
  G4HadronicInteraction* model = mBuilder->GetModel();
  model->SetMinEnergy(emin);
  model->SetMaxEnergy(emax);
  if (verboseLevel>1) {
    G4cout << "### G4VHadronPhysics <" 
	   << model->GetModelName() << " Emin(GeV)= " 
	   << emin/GeV << "  Emax(GeV)= " << emax/GeV
	   << G4endl;
  }

  return model;
}

G4HadronicInteraction* 
G4VHadronPhysics::NewModel(G4HadronicInteraction* model,
			   G4double emin, 
			   G4double emax)
{
  if(!model) return model;
  model->SetMinEnergy(emin);
  model->SetMaxEnergy(emax);
  if (verboseLevel>1) {
    G4cout << "### G4VHadronPhysics <" 
	   << model->GetModelName() << " Emin(GeV)= " 
	   << emin/GeV << "  Emax(GeV)= " << emax/GeV
	   << G4endl;
  }
  return model;
}

