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
// $Id: G4VHadronInelasticPhysics.cc,v 1.1 2009-07-02 09:32:41 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4VHadronInelasticPhysics
//
// Author: 28 June 2009 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4VHadronInelasticPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4HadronicProcessType.hh"
#include "G4HadronicProcessStore.hh"
#include "G4Neutron.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

G4VHadronInelasticPhysics::G4VHadronInelasticPhysics(const G4String& aName, 
						     G4int verb)
  : G4VPhysicsConstructor(aName)
{
  SetVerboseLevel(verb);
  if (verboseLevel>1) {
    G4cout << "### G4VHadronInelasticPhysics: <" << aName << "> is created "
	   << G4endl;
  }
}

G4VHadronInelasticPhysics::~G4VHadronInelasticPhysics() 
{
  G4int n = builders.size();
  if(n > 0) {
    for(G4int i=0; i<n; i++) {delete builders[i];}
  }                           
}                                     

void G4VHadronInelasticPhysics::ConstructParticle()
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

/*
void G4VHadronInelasticPhysics::ConstructProcess()
{
}
*/

G4HadronicInteraction* 
G4VHadronInelasticPhysics::BuildModel(G4VHadronModelBuilder* mBuilder,
				      G4double emin, 
				      G4double emax)
{
  G4HadronicInteraction* model = 0;
  G4int n = builders.size();
  if(n > 0) {
    for(G4int i=0; i<n; i++) {
      if(builders[i] == mBuilder) {
	if (verboseLevel>0) {
	  G4cout << "### G4VHadronInelasticPhysics <" 
		 << mBuilder->GetName() << " already exist "
		 << G4endl;
	}
	model = mBuilder->GetModel();
      }
    }
  }  
  if(!model) {
    builders.push_back(mBuilder);                           
    model =  mBuilder->GetModel();
  }
  model->SetMinEnergy(emin);
  model->SetMaxEnergy(emin);
  if (verboseLevel>1) {
    G4cout << "### G4VHadronInelasticPhysics <" 
	   << model->GetModelName() << " Emin(GeV)= " 
	   << emin/GeV << "  Emax(GeV)= " << emax/GeV
	   << G4endl;
  }

  return model;
}

void 
G4VHadronInelasticPhysics::AddInelasticCrossSection(const G4String& pname, 
						    G4VCrossSectionDataSet* xsec)
{
  const G4ParticleDefinition* p =
    G4ParticleTable::GetParticleTable()->FindParticle(pname);
  if(!p) {
    G4cout << "### G4VHadronInelasticPhysics WARNING: fails to find particle "
	   << pname << G4endl;
  } else {
    AddInelasticCrossSection(p, xsec);
  }
}
 
void 
G4VHadronInelasticPhysics::AddInelasticCrossSection(const G4ParticleDefinition* p, 
						    G4VCrossSectionDataSet* xsec)
{
  G4HadronicProcess* proc =
    G4HadronicProcessStore::Instance()->FindProcess(p, fHadronInelastic);
  if(proc) {
    proc->AddDataSet(xsec);
    if (verboseLevel>1) {
      G4cout << "### G4VHadronInelasticPhysics: the inelastic cross section " 
	     << " is added for " << p->GetParticleName() 
	     << G4endl;
    }
  }
}

void 
G4VHadronInelasticPhysics::AddCaptureCrossSection(G4VCrossSectionDataSet* xsec)
{
  G4HadronicProcess* proc =
    G4HadronicProcessStore::Instance()->FindProcess(G4Neutron::Neutron(), 
						    fCapture);
  if(proc) {
    proc->AddDataSet(xsec);
    if (verboseLevel>1) {
      G4cout << "### G4VHadronInelasticPhysics: the capture cross section " 
	     << " is added for neutron" 
	     << G4endl;
    }
  }
}

void 
G4VHadronInelasticPhysics::AddFissionCrossSection(G4VCrossSectionDataSet* xsec)
{
  G4HadronicProcess* proc =
    G4HadronicProcessStore::Instance()->FindProcess(G4Neutron::Neutron(), 
						    fFission);
  if(proc) {
    proc->AddDataSet(xsec);
    if (verboseLevel>1) {
      G4cout << "### G4VHadronInelasticPhysics: the capture cross section " 
	     << " is added for neutron" 
	     << G4endl;
    }
  }
}

