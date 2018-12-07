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
/// \file hadronic/Hadr02/src/IonHIJINGPhysics.cc
/// \brief Implementation of the IonHIJINGPhysics class
//
//
//---------------------------------------------------------------------------
//
// Class:    IonHIJINGPhysics
//
// Author:   2012 Andrea Dotti
//
//
// Modified:
//
// ------------------------------------------------------------
// 
#ifdef G4_USE_HIJING
#include "IonHIJINGPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronInelasticProcess.hh"

#include "G4HIJING_Model.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4BuilderType.hh"
#include "G4HadronicParameters.hh"
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonHIJINGPhysics::IonHIJINGPhysics(G4int ver)
  : G4VHadronPhysics("ionInelasticHIJING"),fVerbose(ver),
    fWasActivated(false)
{
  fModel = 0;
  SetPhysicsType(bIons);
  if(fVerbose > 1) { G4cout << "### IonHIJINGPhysics" << G4endl; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonHIJINGPhysics::~IonHIJINGPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void IonHIJINGPhysics::ConstructProcess()
{
  if(fWasActivated) { return; }
  fWasActivated = true;

  G4double emin = 0.*MeV;
  G4double emaxFTF = 25.*GeV;
  G4double eminHIJ = 12.*GeV;
  G4double emaxHIJ = G4HadronicParameters::Instance()->GetMaxEnergy();

  G4ExcitationHandler* handler = new G4ExcitationHandler();
  G4PreCompoundModel* thePreCompound = new G4PreCompoundModel(handler);

  // Binary Cascade
  theIonBC = new G4BinaryLightIonReaction(thePreCompound);
  theIonBC->SetMinEnergy(0.0);
  theIonBC->SetMaxEnergy(4*GeV); //4

  // FTFP
  theBuilder = new G4FTFBuilder("FTFP",thePreCompound);
  theFTFP = theBuilder->GetModel();
  theFTFP->SetMinEnergy(2*GeV);
  theFTFP->SetMaxEnergy(emaxFTF);
  
  //HIJING
  fModel = new G4HIJING_Model();
  fModel->SetMinEnergy( eminHIJ );
  fModel->SetMaxEnergy( emaxHIJ );

  theNuclNuclData = 
    new G4CrossSectionInelastic( new G4ComponentGGNuclNuclXsc() );

  AddProcess("dInelastic", G4Deuteron::Deuteron(),false);
  AddProcess("tInelastic",G4Triton::Triton(),false);
  AddProcess("He3Inelastic",G4He3::He3(),true);
  AddProcess("alphaInelastic", G4Alpha::Alpha(),true);
  AddProcess("ionInelastic",G4GenericIon::GenericIon(),true);

  if(fVerbose > 1) {
    G4cout << "IonHIJINGPhysics::ConstructProcess done! " 
   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IonHIJINGPhysics::AddProcess(const G4String& name, 
                                  G4ParticleDefinition* part,  G4bool isIon)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, part);
  G4ProcessManager* pManager = part->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);
  hadi->AddDataSet(theNuclNuclData);
  hadi->RegisterMe( theIonBC );
  hadi->RegisterMe( theFTFP );
  hadi->RegisterMe( fModel );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif //HIJING
