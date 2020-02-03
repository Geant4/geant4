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
// Class:    G4IonPhysicsXS
//
// Author:      A.Ivanchenko 28.07.2018
//
// Modified: 
//
//---------------------------------------------------------------------------

#include "G4IonPhysicsXS.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4GenericIon.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4ParticleInelasticXS.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4HadronicInteraction.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"

using namespace std;

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4IonPhysicsXS);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonPhysicsXS::G4IonPhysicsXS(G4int ver)
  : G4IonPhysicsXS("ionPhysicsXS")
{
  verbose = ver;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonPhysicsXS::G4IonPhysicsXS(const G4String& nname)
  : G4IonPhysics(nname)
{
  G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
  param->SetDeexChannelsType(fCombined);
  if(verbose > 1) { G4cout << "### IonPhysics: " << nname << G4endl; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IonPhysicsXS::~G4IonPhysicsXS()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4IonPhysicsXS::AddProcess(const G4String& name, 
		 	        G4ParticleDefinition* part,
			        G4HadronicInteraction* theIonBC,
			        G4HadronicInteraction* theFTFP,
			        G4VCrossSectionDataSet* xs)
{
  G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, part);
  G4ProcessManager* pManager = part->GetProcessManager();
  pManager->AddDiscreteProcess(hadi);
  
  if(part == G4GenericIon::GenericIon()) {  
    hadi->AddDataSet(xs);
  } else {
    hadi->AddDataSet(new G4ParticleInelasticXS(part));
  }  
  hadi->RegisterMe(theIonBC);
  hadi->RegisterMe(theFTFP);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
