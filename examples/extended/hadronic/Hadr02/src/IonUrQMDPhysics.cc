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
/// \file hadronic/Hadr02/src/IonUrQMDPhysics.cc
/// \brief Implementation of the IonUrQMDPhysics class
//
//
//---------------------------------------------------------------------------
//
// Class:    IonUrQMDPhysics
//
// Author:   2012 Andrea Dotti
//
//
// Modified:
//
// ------------------------------------------------------------
// 
#ifdef G4_USE_URQMD
#include "IonUrQMDPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4TripathiCrossSection.hh"
#include "G4TripathiLightCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4IonProtonCrossSection.hh"

#include "G4UrQMD1_3Model.hh"

#include "G4BuilderType.hh"
#include "G4HadronicParameters.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonUrQMDPhysics::IonUrQMDPhysics(G4int ver)
  : G4VHadronPhysics("ionInelasticUrQMD"),verbose(ver),
    fWasActivated(false)
{
  fTripathi = fTripathiLight = fShen = fIonH = 0;
  fModel = 0;
  SetPhysicsType(bIons);
  if(fVerbose > 1) { G4cout << "### IonUrQMDPhysics" << G4endl; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IonUrQMDPhysics::~IonUrQMDPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void IonUrQMDPhysics::ConstructProcess()
{
  if(fWasActivated) { return; }
  fWasActivated = true;

  G4double emin = 0.*MeV;
  G4double emax = G4HadronicParameters::Instance()->GetMaxEnergy();

  fModel = new G4UrQMD1_3Model();
  fModel->SetMinEnergy( emin );
  fModel->SetMaxEnergy( emax );

  fShen = new G4IonsShenCrossSection();
  fTripathi = new G4TripathiCrossSection();
  fTripathiLight = new G4TripathiLightCrossSection();
  fIonH = new G4IonProtonCrossSection();
  fShen->SetMaxKinEnergy( emax );
  fTripathi->SetMaxKinEnergy( emax );
  fTripathiLight->SetMaxKinEnergy( emax );
  fIonH->SetMaxKinEnergy( emax );


  AddProcess("dInelastic", G4Deuteron::Deuteron(),false);
  AddProcess("tInelastic",G4Triton::Triton(),false);
  AddProcess("He3Inelastic",G4He3::He3(),true);
  AddProcess("alphaInelastic", G4Alpha::Alpha(),true);
  AddProcess("ionInelastic",G4GenericIon::GenericIon(),true);

  if(fVerbose > 1) {
    G4cout << "IonUrQMDPhysics::ConstructProcess done! " 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IonUrQMDPhysics::AddProcess(const G4String& name, 
                                 G4ParticleDefinition* part, 
                                 G4bool isIon)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, part);
  G4ProcessManager* pManager = part->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);
  hadi->AddDataSet(fShen);
  //  hadi->AddDataSet(fTripathi);
  // hadi->AddDataSet(fTripathiLight);
  if(isIon) { hadi->AddDataSet(fIonH); }
  hadi->RegisterMe( fModel );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif //URQMD
