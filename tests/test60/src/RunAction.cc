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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction() : fpProcCounter(0)
{  
	fpPrimary = 0;
	fTotalCount = 0;
	fE_Transfered = 0.;
	fInitialized = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{
	delete fpProcCounter;
}

void RunAction::Initialize()
{
	fpDetector = (DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
	fpPrimary = (PrimaryGeneratorAction*) G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
	fpProcCounter = new ProcessesCount;
	fInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run*)
{  
	if(fInitialized == false) Initialize();
	else fpProcCounter->clear();

	fE_Transfered = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcesses(G4String procName)
{
	//does the process  already encounted ?
	size_t nbProc = fpProcCounter->size();
	size_t i = 0;
	while ((i<nbProc)&&((*fpProcCounter)[i]->GetName()!=procName)) i++;
	if (i == nbProc) fpProcCounter->push_back( new OneProcessCount(procName));

	(*fpProcCounter)[i]->Count();
}

void RunAction::ResetCounter()
{
	// delete and remove all contents in ProcCounter
	while (fpProcCounter->size()>0){
		OneProcessCount* aProcCount=fpProcCounter->back();
		fpProcCounter->pop_back();
		delete aProcCount;
	}
	// delete ProcCounter;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run* aRun)
{
	G4int NbOfEvents = aRun->GetNumberOfEvent();
	if (NbOfEvents == 0) return;

	G4int  prec = G4cout.precision(5);

	G4Material* material = fpDetector->GetMaterial();
	G4double density = material->GetDensity();
	G4int survive = 0;

	const G4ParticleDefinition* particle = fpPrimary->GetParticleGun()->GetParticleDefinition();
	const G4String& particleName = particle->GetParticleName();
	G4double energy = fpPrimary->GetParticleGun()->GetParticleEnergy();

	G4cout << "\n The run consists of " << NbOfEvents << " "<< particleName << " of "
			<< G4BestUnit(energy,"Energy") << " through "
			<< G4BestUnit(fpDetector->GetSize(),"Length") << " of "
			<< material->GetName() << " (density: "
			<< G4BestUnit(density,"Volumic Mass") << ")" << G4endl;

	//frequency of processes
	G4cout << "\n Process calls frequency --->";

	for (size_t i=0; i< fpProcCounter->size();i++)
	{
		G4String procName = (*fpProcCounter)[i]->GetName();
		G4int    count    = (*fpProcCounter)[i]->GetCounter();
		G4cout << "\t" << procName << " = " << count;
		if (procName == "Transportation") survive = count;
	}

	if (survive > 0) {
		G4cout << "\n\n Nb of incident particles surviving after "
				<< G4BestUnit(fpDetector->GetSize(),"Length") << " of "
				<< material->GetName() << " : " << survive << G4endl;
	}

	if (fTotalCount == 0) fTotalCount = 1;   //force printing anyway

	G4cout << G4endl;
	G4cout << " Total transfered energy (keV)=" <<fE_Transfered/keV<< G4endl;
	G4cout << G4endl;

	/*
  //compute mean free path and related quantities
  //
  G4double MeanFreePath = sumTrack /totalCount;
  G4double MeanTrack2   = sumTrack2/totalCount;
  G4double rms = std::sqrt(std::fabs(MeanTrack2 - MeanFreePath*MeanFreePath));
  G4double CrossSection = 1./MeanFreePath;
  G4double massicMFP = MeanFreePath*density;
  G4double massicCS  = 1./massicMFP;

  G4cout << "\n\n MeanFreePath:\t"   << G4BestUnit(MeanFreePath,"Length")
         << " +- "                   << G4BestUnit( rms,"Length")
         << "\tmassic: "             << G4BestUnit(massicMFP, "Mass/Surface")
         << "\n CrossSection:\t"     << CrossSection*cm << " cm^-1 "
	 << "\t\t\tmassic: "         << G4BestUnit(massicCS, "Surface/Mass")
         << G4endl;

  //compute energy transfer coefficient
  //
  G4double MeanTransfer   = eTransfer/totalCount;
  G4double massTransfCoef = massicCS*MeanTransfer/energy;

  G4cout << "\n mean energy of charged secondaries: " << G4BestUnit(MeanTransfer, "Energy")
	 << "\tmass_energy_transfer coef: "           << G4BestUnit(massTransfCoef, "Surface/Mass")
         << G4endl;

  //check cross section from G4EmCalculator
  //
  G4cout << "\n Verification : "
         << "crossSections from G4EmCalculator \n";

  G4EmCalculator emCalculator;
  G4double sumc = 0.0;
  for (size_t i=0; i< ProcCounter->size();i++) {
    G4String procName = (*ProcCounter)[i]->GetName();
    G4double massSigma =
    emCalculator.GetCrossSectionPerVolume(energy,particle,
                                              procName,material)/density;
    if (particle == G4Gamma::Gamma())
       massSigma =
       emCalculator.ComputeCrossSectionPerVolume(energy,particle,
                                              procName,material)/density;
    sumc += massSigma;
    G4cout << "\t" << procName << "= "
           << G4BestUnit(massSigma, "Surface/Mass");
  }
  G4cout << "\ttotal= "
         << G4BestUnit(sumc, "Surface/Mass") << G4endl;

	 */
	//restore default format
	G4cout.precision(prec);
	ResetCounter();

}
