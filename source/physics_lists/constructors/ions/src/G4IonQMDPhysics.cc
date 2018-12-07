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
// ClassName:   G4IonQMDPhysics
//   Created from G4IonBinaryCascadePhysics
//
// Author:      G.Folger
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4IonQMDPhysics.hh"

#include "G4SystemOfUnits.hh"

#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4IonConstructor.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4QMDReaction.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4BuilderType.hh"

#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Nuclei
#include "G4IonConstructor.hh"
#include "G4BuilderType.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicParameters.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonQMDPhysics);

G4ThreadLocal G4FTFBuilder* G4IonQMDPhysics::theFTFPBuilder = nullptr;

G4IonQMDPhysics::G4IonQMDPhysics(G4int ver)
  :  G4IonQMDPhysics("IonQMD", ver)
{}

G4IonQMDPhysics::G4IonQMDPhysics(const G4String& nname, G4int ver)
  :  G4VPhysicsConstructor(nname), verbose(ver)
{
  eminQMD  = 100.*MeV;
  emaxQMD  = 10.*GeV;
  overlap  = 10*MeV;
  SetPhysicsType(bIons);
  G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
  param->SetDeexChannelsType(fCombined);
  if(verbose > 1) { G4cout << "### IonPhysics: " << nname << G4endl; }
}

G4IonQMDPhysics::~G4IonQMDPhysics()
{
  delete theFTFPBuilder; theFTFPBuilder = nullptr;
}

void G4IonQMDPhysics::ConstructProcess()
{
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  G4PreCompoundModel* thePreCompound = static_cast<G4PreCompoundModel*>(p); 
  if(!thePreCompound) { thePreCompound = new G4PreCompoundModel; }

  G4BinaryLightIonReaction* theIonBC = new G4BinaryLightIonReaction(thePreCompound);
  theIonBC->SetMaxEnergy(eminQMD + overlap);

  G4double emax = G4HadronicParameters::Instance()->GetMaxEnergy();
  G4HadronicInteraction* theFTFP = nullptr;
  if(emax > emaxQMD) {
    theFTFPBuilder = new G4FTFBuilder("FTFP",thePreCompound);
    theFTFP = theFTFPBuilder->GetModel();
    theFTFP->SetMinEnergy(emaxQMD - overlap);
    theFTFP->SetMaxEnergy(emax);
  }

  G4QMDReaction* theQMD = new G4QMDReaction();
  theQMD->SetMinEnergy(eminQMD);
  theQMD->SetMaxEnergy(emaxQMD);

  G4VCrossSectionDataSet* theNuclNuclData = 
    new G4CrossSectionInelastic( new G4ComponentGGNuclNuclXsc() );

  AddProcess("dInelastic", G4Deuteron::Deuteron(), theIonBC, theQMD, theFTFP, theNuclNuclData);
  AddProcess("tInelastic", G4Triton::Triton(), theIonBC, theQMD, theFTFP, theNuclNuclData);
  AddProcess("He3Inelastic", G4He3::He3(), theIonBC, theQMD, theFTFP, theNuclNuclData);
  AddProcess("alphaInelastic", G4Alpha::Alpha(), theIonBC, theQMD, theFTFP, theNuclNuclData);
  AddProcess("ionInelastic", G4GenericIon::GenericIon(), theIonBC, theQMD, theFTFP, theNuclNuclData);
}

void G4IonQMDPhysics::AddProcess(const G4String& name,
			         G4ParticleDefinition* p, 
				 G4BinaryLightIonReaction* BIC,
				 G4QMDReaction* QMD,
				 G4HadronicInteraction* FTFP,
				 G4VCrossSectionDataSet* theNuclNuclData)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, p);
  G4ProcessManager* pManager = p->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);

  hadi->AddDataSet(theNuclNuclData);

  hadi->RegisterMe(BIC);
  hadi->RegisterMe(QMD);
  if(FTFP) { hadi->RegisterMe(FTFP); }  

  if(verbose > 1) {
    G4cout << "Register " << hadi->GetProcessName()
	   << " for " << p->GetParticleName() << G4endl
	   << "       Binary Cascade for E(MeV)= 0 - " 
	   << eminQMD+overlap;
    G4cout << "       QMD for E(MeV)= " << eminQMD << " - " << emaxQMD;
    if(FTFP) {
      G4cout << "       FTFP for E(MeV)= " << emaxQMD-overlap << " - " << FTFP->GetMaxEnergy();
    }
    G4cout << G4endl;
  }
}

void G4IonQMDPhysics::ConstructParticle()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}
