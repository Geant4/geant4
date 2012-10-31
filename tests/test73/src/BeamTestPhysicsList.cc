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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "BeamTestPhysicsList.hh"

#include "G4EmStandardPhysics_option1LHCb_M.hh" // LHCb below code plus additions by MR ionisation and step size 
#include "MyEmPhysicsList.hh" // Private study
#include "PhysListEmStandard.hh" // LHCb code
#include "G4SystemOfUnits.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BeamTestPhysicsList::BeamTestPhysicsList() 
: G4VModularPhysicsList()
{
	G4LossTableManager::Instance();

	currentDefaultCut   = 1.0*mm;
	cutForGamma         = currentDefaultCut;
	cutForElectron      = currentDefaultCut;
	cutForPositron      = currentDefaultCut;
	cutForMuon     	    = currentDefaultCut;
	cutForPion    	    = currentDefaultCut;


	SetVerboseLevel(1);
	// EM physics
	//name = "emstandard_opt0";
	//name = "emstandard_opt1";
	//name = "emstandard_opt2";
	//name = "emstandard_opt3";
	//name = "empenelope";
	//name = "emlivermore";
	name = G4String("local");
	//name = G4String("LHCb");
	//name = G4String("LHCb_M");
	//emPhysicsList = new PhysListEmStandard();
	emPhysicsList = new MyEmPhysicsList(name);/*G4EmStandardPhysics_option2();*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BeamTestPhysicsList::~BeamTestPhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

// leptons
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

// Mesons
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Eta.hh"
#include "G4EtaPrime.hh"

#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4AntiKaonZero.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"

// Baryons
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"

// Nuclei
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamTestPhysicsList::ConstructParticle()
{
	// pseudo-particles
	G4Geantino::GeantinoDefinition();
	G4ChargedGeantino::ChargedGeantinoDefinition();

	// gamma
	G4Gamma::GammaDefinition();

	// optical photon
	G4OpticalPhoton::OpticalPhotonDefinition();

	// leptons
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
	G4MuonPlus::MuonPlusDefinition();
	G4MuonMinus::MuonMinusDefinition();

	G4NeutrinoE::NeutrinoEDefinition();
	G4AntiNeutrinoE::AntiNeutrinoEDefinition();
	G4NeutrinoMu::NeutrinoMuDefinition();
	G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();  

	// mesons
	G4PionPlus::PionPlusDefinition();
	G4PionMinus::PionMinusDefinition();
	G4PionZero::PionZeroDefinition();
	G4Eta::EtaDefinition();
	G4EtaPrime::EtaPrimeDefinition();
	G4KaonPlus::KaonPlusDefinition();
	G4KaonMinus::KaonMinusDefinition();
	G4KaonZero::KaonZeroDefinition();
	G4AntiKaonZero::AntiKaonZeroDefinition();
	G4KaonZeroLong::KaonZeroLongDefinition();
	G4KaonZeroShort::KaonZeroShortDefinition();

	// barions
	G4Proton::ProtonDefinition();
	G4AntiProton::AntiProtonDefinition();
	G4Neutron::NeutronDefinition();
	G4AntiNeutron::AntiNeutronDefinition();

	// ions
	G4Deuteron::DeuteronDefinition();
	G4Triton::TritonDefinition();
	G4Alpha::AlphaDefinition();
	G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4EmProcessOptions.hh"

void BeamTestPhysicsList::ConstructProcess()
{
	// Transportation
	//
	AddTransportation();

	// Electromagnetic physics list
	//
	AddPhysicsList();
	emPhysicsList->ConstructProcess();
	/*	if (name == "local") 
		{
		emPhysicsList->AddStepMax();
		}*/

	// Em options
	//
	// Main options and setting parameters are shown here.
	// Several of them have default values.
	//
	G4EmProcessOptions emOptions;

	//physics tables
	//
	//emOptions.SetMinEnergy(100*eV);	//default    
	//emOptions.SetMaxEnergy(100*TeV);	//default  
	//emOptions.SetDEDXBinning(12*20);	//default=12*7  
	//emOptions.SetLambdaBinning(12*20);	//default=12*7

	emOptions.SetBuildCSDARange(true);     
	//emOptions.SetMaxEnergyForCSDARange(100*TeV);
	//emOptions.SetDEDXBinningForCSDARange(12*20);

	//emOptions.SetSplineFlag(true);	//default

	emOptions.SetVerbose(0);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamTestPhysicsList::AddPhysicsList()
{
	//G4String name;
	//if (verboseLevel>0) {
	G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
	//}

	//if (name == 0) return;

	if (name == "LHCb") {

		delete emPhysicsList;
		emPhysicsList = new PhysListEmStandard();

	} else if (name == "LHCb_M"){
		delete emPhysicsList;
		emPhysicsList = new G4EmStandardPhysics_option1LHCb_M();

	}else if (name == "local"){
		delete emPhysicsList;
		emPhysicsList = new MyEmPhysicsList(name);

	}else if (name == "emstandard_opt0"){
		delete emPhysicsList;
		emPhysicsList = new G4EmStandardPhysics();

	} else if (name == "emstandard_opt1"){
		delete emPhysicsList;
		emPhysicsList = new G4EmStandardPhysics_option1();

	} else if (name == "emstandard_opt2"){
		delete emPhysicsList;
		emPhysicsList = new G4EmStandardPhysics_option2();

	} else if (name == "emstandard_opt3"){
		delete emPhysicsList;
		emPhysicsList = new G4EmStandardPhysics_option3();

	} else if (name == "empenelope"){
		delete emPhysicsList;
		emPhysicsList = new G4EmPenelopePhysics();

	} else if (name == "emlivermore"){
		delete emPhysicsList;
		emPhysicsList = new G4EmLivermorePhysics();

	} else {

		G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
			<< " is not defined"
			<< G4endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

void BeamTestPhysicsList::SetCuts()
{ 
	// fixe lower limit for cut
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);

	// set cut values for gamma at first and for e- second and next for e+,
	// because some processes for e+/e- need cut values for gamma
	SetCutValue(cutForGamma, "gamma");
	SetCutValue(cutForElectron, "e-");
	SetCutValue(cutForPositron, "e+");
	//SetCutValue(cutForMuon, "mu-");
	//SetCutValue(cutForPion, "pi-");
	DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void BeamTestPhysicsList::SetCutForGamma(G4double cut)
  {
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamTestPhysicsList::SetCutForElectron(G4double cut)
{
cutForElectron = cut;
SetParticleCuts(cutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamTestPhysicsList::SetCutForPositron(G4double cut)
{
cutForPositron = cut;
SetParticleCuts(cutForPositron, G4Positron::Positron());
}
 */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


