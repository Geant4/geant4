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
// $Id: G4VHadronPhysics.cc 85579 2014-10-31 09:04:00Z gcosmo $
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
#include "G4VCrossSectionDataSet.hh"
#include "G4HadronicProcessType.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4Neutron.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"

G4ThreadLocal std::vector<G4VHadronModelBuilder*>* G4VHadronPhysics::builders = 0;

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
{
    if ( builders ) {
        G4int n = builders->size();
        if(n > 0) {
            for(G4int i=0; i<n; i++) {delete (*builders)[i];}
        }
        delete builders;
    }
}

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
  if ( builders == 0 ) builders = new std::vector<G4VHadronModelBuilder*>;
  builders->push_back(mBuilder);                           
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

void 
G4VHadronPhysics::AddInelasticCrossSection(const G4String& pname, 
					   G4VCrossSectionDataSet* xsec)
{
  const G4ParticleDefinition* p =
    G4ParticleTable::GetParticleTable()->FindParticle(pname);
  if(!p) {
    G4cout << "### G4VHadronPhysics WARNING: fails to find particle "
	   << pname << G4endl;
  } else {
    AddInelasticCrossSection(p, xsec);
  }
}
 
void 
G4VHadronPhysics::AddInelasticCrossSection(const G4ParticleDefinition* p, 
					   G4VCrossSectionDataSet* xsec)
{
  if(!p) return;
  G4HadronicProcess* had = FindInelasticProcess(p);
  if(!had) return;
  had->AddDataSet(xsec);
  if (verboseLevel>1) {
    G4cout << "### G4VHadronPhysics: the inelastic cross section " 
	   << " is added for " << p->GetParticleName() 
	   << G4endl;
  }
}

void 
G4VHadronPhysics::AddElasticCrossSection(const G4String& pname, 
					   G4VCrossSectionDataSet* xsec)
{
  const G4ParticleDefinition* p =
    G4ParticleTable::GetParticleTable()->FindParticle(pname);
  if(!p) {
    G4cout << "### G4VHadronPhysics WARNING: fails to find particle "
	   << pname << G4endl;
  } else {
    AddElasticCrossSection(p, xsec);
  }
}
 
void 
G4VHadronPhysics::AddElasticCrossSection(const G4ParticleDefinition* p, 
					   G4VCrossSectionDataSet* xsec)
{
  if(!p) return;
  G4HadronicProcess* had = FindElasticProcess(p);
  if(!had) return;
  had->AddDataSet(xsec);
  if (verboseLevel>1) {
    G4cout << "### G4VHadronPhysics: the inelastic cross section " 
	   << " is added for " << p->GetParticleName() 
	   << G4endl;
  }
}

void 
G4VHadronPhysics::AddCaptureCrossSection(G4VCrossSectionDataSet* xsec)
{
  G4HadronicProcess* had = FindCaptureProcess();
  if(!had) return;
  had->AddDataSet(xsec);
  if (verboseLevel>1) {
    G4cout << "### G4VHadronPhysics: the capture cross section " 
	   << " is added for neutron" 
	   << G4endl;
  }
}

void 
G4VHadronPhysics::AddFissionCrossSection(G4VCrossSectionDataSet* xsec)
{
  G4HadronicProcess* had = FindFissionProcess();
  if(!had) return;
  had->AddDataSet(xsec);
  if (verboseLevel>1) {
    G4cout << "### G4VHadronPhysics: the fission cross section " 
	   << " is added for neutron" 
	   << G4endl;
  }
}

G4HadronicProcess* 
G4VHadronPhysics::FindInelasticProcess(const G4String& pname)
{
  G4HadronicProcess* had = 0;
  const G4ParticleDefinition* p =
    G4ParticleTable::GetParticleTable()->FindParticle(pname);
  if(!p) {
    G4cout << "### G4VHadronPhysics WARNING: fails to find particle "
	   << pname << G4endl;
    return had;
  }
  return FindInelasticProcess(p);
}

G4HadronicProcess* 
G4VHadronPhysics::FindInelasticProcess(const G4ParticleDefinition* p)
{
  G4HadronicProcess* had = 0;
  if(!p) return had;
  G4ProcessManager* pmanager = p->GetProcessManager();
  G4ProcessVector*  pv = pmanager->GetProcessList();
  size_t n = pv->size();
  if(0 < n) {
    for(size_t i=0; i<n; ++i) {
      if(fHadronInelastic == ((*pv)[i])->GetProcessSubType()) {
	had = static_cast<G4HadronicProcess*>((*pv)[i]);
	return had;
      }
    }
  }
  G4ParticleDefinition* part = const_cast<G4ParticleDefinition*>(p);
  had = new G4HadronInelasticProcess(part->GetParticleName()+"Inelastic",part);
  pmanager->AddDiscreteProcess(had);
  return had;
}

G4HadronicProcess* 
G4VHadronPhysics::FindElasticProcess(const G4String& pname)
{
  G4HadronicProcess* had = 0;
  const G4ParticleDefinition* p =
    G4ParticleTable::GetParticleTable()->FindParticle(pname);
  if(!p) {
    G4cout << "### G4VHadronPhysics WARNING: fails to find particle "
	   << pname << G4endl;
    return had;
  }
  return FindElasticProcess(p);
}

G4HadronicProcess* 
G4VHadronPhysics::FindElasticProcess(const G4ParticleDefinition* p)
{
  G4HadronicProcess* had = 0;
  if(!p) return had;
  G4ProcessManager* pmanager = p->GetProcessManager();
  G4ProcessVector*  pv = pmanager->GetProcessList();
  size_t n = pv->size();
  if(0 < n) {
    for(size_t i=0; i<n; ++i) {
      if(fHadronElastic == ((*pv)[i])->GetProcessSubType()) {
	had = static_cast<G4HadronicProcess*>((*pv)[i]);
	return had;
      }
    }
  }
  had = new G4HadronElasticProcess("hElastic");
  pmanager->AddDiscreteProcess(had);
  return had;
}

G4HadronicProcess* G4VHadronPhysics::FindCaptureProcess()
{
  G4HadronicProcess* had = 0;
  G4ProcessManager* pmanager = 
    G4Neutron::Neutron()->GetProcessManager();
  G4ProcessVector*  pv = pmanager->GetProcessList();
  size_t n = pv->size();
  if(0 < n) {
    for(size_t i=0; i<n; ++i) {
      if(fCapture == ((*pv)[i])->GetProcessSubType()) {
	had = static_cast<G4HadronicProcess*>((*pv)[i]);
	return had;
      }
    }
  }
  had = new G4HadronCaptureProcess("nCapture");
  pmanager->AddDiscreteProcess(had);
  return had;
}

G4HadronicProcess* G4VHadronPhysics::FindFissionProcess()
{
  G4HadronicProcess* had = 0;
  G4ProcessManager* pmanager = 
    G4Neutron::Neutron()->GetProcessManager();
  G4ProcessVector*  pv = pmanager->GetProcessList();
  size_t n = pv->size();
  if(0 < n) {
    for(size_t i=0; i<n; ++i) {
      if(fFission == ((*pv)[i])->GetProcessSubType()) {
	had = static_cast<G4HadronicProcess*>((*pv)[i]);
	return had;
      }
    }
  }
  had = new G4HadronFissionProcess("nFission");
  pmanager->AddDiscreteProcess(had);
  return had;
}

