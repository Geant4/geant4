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
// RunAction.cc
// 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, 
PhysicsList* phy, PrimaryGeneratorAction* kin)
:detector(det), physics(phy)  , primary(kin)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //initialisation
  Nsteps=0;
  Nsteps2=0;

  theta=0;
  theta2=0; 
  N_rec=0;
  EnergyDeposit  = EnergyDeposit2  = 0.;
  NonIonEnergyDeposit  = NonIonEnergyDeposit2  = 0.;
  sum_TL=sum_TL2=0.;
  sum_T=sum_T2=0.;

  TrakLenPrim = TrakLenPrim2 = 0.;

  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include <iostream>

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  // compute mean and rms
  //
  G4int TotNbofEvents = aRun->GetNumberOfEvent();
  //G4cout << "Ntot= " << TotNbofEvents << G4endl;
  if (TotNbofEvents == 0) return;
  
  //total energy loss
  EnergyDeposit /= TotNbofEvents; EnergyDeposit2 /= TotNbofEvents;
  G4double rmsEdep = EnergyDeposit2 - EnergyDeposit*EnergyDeposit;
  if(rmsEdep >0.)rmsEdep=std::sqrt(rmsEdep/TotNbofEvents);
  else rmsEdep =0;

  //nuclear energy loss
  NonIonEnergyDeposit /= TotNbofEvents; NonIonEnergyDeposit2 /= TotNbofEvents;
  G4double rmsEnondep = NonIonEnergyDeposit2 - NonIonEnergyDeposit*NonIonEnergyDeposit;
  if(rmsEnondep>0.) rmsEnondep= std::sqrt(rmsEnondep/TotNbofEvents);
  else rmsEnondep=0;

  //mean sum of T( kinetic energy of secondary)x L(T) partition energy 
  sum_TL/=TotNbofEvents;     sum_TL2/=TotNbofEvents;
  G4double rmssum_TL =sum_TL2- sum_TL*sum_TL;
  if(rmssum_TL>0.) rmssum_TL=std::sqrt(rmssum_TL/TotNbofEvents);
  else rmssum_TL =0;

  //mean kinetic energy of secondary particles (IDp==1) 
  G4double rmssum_T = 0.0;
  if(N_rec > 0) {
    sum_T/=N_rec;     sum_T2/=N_rec;
    G4double rmssum_T =sum_T2- sum_T*sum_T;
    if(rmssum_T>0.) rmssum_T=std::sqrt(rmssum_T/N_rec);
  }

  //mean number of steps:
  Nsteps/=TotNbofEvents;  Nsteps2/=TotNbofEvents;
  G4double rmsSteps= Nsteps2 -Nsteps*Nsteps;
  if(rmsSteps>0) rmsSteps= std::sqrt(rmsSteps/TotNbofEvents);
  else rmsSteps=0;

  //scattering angle
  theta/=TotNbofEvents ; theta2/=TotNbofEvents;	
  G4double rmsTheta =theta2-theta*theta;
  if(rmsTheta>0.) rmsTheta =std::sqrt(rmsTheta/TotNbofEvents);
  else rmsTheta =0;
  
  //track length
  TrakLenPrim /= TotNbofEvents; TrakLenPrim2 /= TotNbofEvents;
  G4double rmsTLPrim = TrakLenPrim2 - TrakLenPrim*TrakLenPrim;
  if (rmsTLPrim>0.) rmsTLPrim = std::sqrt(rmsTLPrim/TotNbofEvents);
  else rmsTLPrim = 0.;

  //..............................................................

  G4Material* material = detector->GetAbsorberMaterial();
  G4double thickness  = detector->GetAbsorberThickness();
  G4double density = material->GetDensity();
   
  G4ParticleDefinition* particle = primary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double energy = primary->GetParticleGun()->GetParticleEnergy();
   
  //Stopping Power and NIEL from simulation.
  //
  // effective length
  G4double length=TrakLenPrim;

  //G4cout << "Length= " << length << G4endl;
  //G4cout << "Density= " << density << G4endl;
  //G4cout << "Nsteps= " << Nsteps << G4endl;

  // total energy loss  
  G4double meandEdx  = EnergyDeposit/length;
  // nuclear energy loss
  G4double meandEdx_nucl  = NonIonEnergyDeposit/length;
  // NIEL 
  G4double meandEdx_sumTL=sum_TL/length;

	//[MeVcm2/g]
  G4double stopPower = meandEdx/density;  
  G4double stopPower_nucl = meandEdx_nucl/density;
  G4double stopPower_sumTL=meandEdx_sumTL/density;

 
 	//mean free path & corss section 
  G4double freepath= TrakLenPrim/Nsteps;
  G4double er1=rmsTLPrim/Nsteps;	
  G4double er2=freepath*rmsSteps/Nsteps;
  G4double rmsFreepath=std::sqrt(er1*er1+er2*er2);


  G4double NA  = material->GetTotNbOfAtomsPerVolume();
  //G4cout << "NA= " << NA << G4endl;
  G4double CrossSection =1./(NA*freepath); 
  G4double rmsCrossSection=rmsFreepath*CrossSection/freepath;

  Th=physics->GetEth();

  G4cout << "\n ======================== run summary ======================\n";

  G4cout.precision(4);
  
  G4cout << "\n The run was " << TotNbofEvents << " " << partName << " of "
         << G4BestUnit(energy,"Energy") 
         << ", through " 
	 << G4BestUnit(thickness,"Length") << " of "
	 << material->GetName() << "\n (density: " 
	 << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  G4cout << "\n Threshold Energy for displacement \n "
        << Th/eV<<" eV" <<G4endl;


  G4cout << "\n ============= Primary particle statistics ==============\n";
 
  G4cout << "\n Total track length in absorber:\n "
         << G4BestUnit(TrakLenPrim,"Length") << " +- "
         << G4BestUnit(rmsTLPrim,       "Length") << G4endl;


  G4cout << "\n Mean Number of Steps :\n "<<Nsteps<< " +/- "<< rmsSteps<<G4endl;
								
  G4cout << "\n Mean Free Path :\n "<<G4BestUnit(freepath,"Length")<<  
				" +/- "<< G4BestUnit(rmsFreepath,"Length")<<G4endl;
  
  G4cout << "\n Mean Cross Section :\n "<<G4BestUnit(CrossSection,"Surface")<<
			" +/- "<<G4BestUnit(rmsCrossSection,"Surface")<<G4endl;

  G4cout << "\n Mean scattering angle :\n "<<G4BestUnit(theta,"Angle")<< " +/- "
			<< G4BestUnit(rmsTheta,"Angle")<<G4endl;

  
  G4cout << "\n Total energy deposit in absorber:\n "
         << G4BestUnit(EnergyDeposit,"Energy") << " +/- "
         << G4BestUnit(rmsEdep,      "Energy") 
         << G4endl;
  G4cout << "-----> dE/dx total= " << meandEdx/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;


  G4cout << "\n Nuclear energy deposit in absorber:\n "
         << G4BestUnit(NonIonEnergyDeposit,"Energy") << " +/- "
         << G4BestUnit(rmsEnondep,      "Energy")
         << G4endl;
  G4cout << "-----> dE/dx  nucl = " << meandEdx_nucl/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower_nucl/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;


  G4cout <<"\n NIEL in absorber (Th>"<<Th/eV <<" eV):\n "
         << G4BestUnit(sum_TL,"Energy") << " +/- "
         << G4BestUnit(rmssum_TL,      "Energy")
         << G4endl;
  G4cout << "-----> NIEL = " << meandEdx_sumTL/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower_sumTL/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;

  //Frenkel-pairs Concentration:
  G4double l=	detector->GetAbsorberSizeYZ();
  G4double area=l*l;//cm2
  G4double flu =TotNbofEvents/area; // #/cm2
  G4double Edis= meandEdx_sumTL*flu;// MeV/cm3
  G4double FP= Edis/(2.5*Th);
 
  G4cout<<"\n----------------------------------------------------------------"<<G4endl;
  G4cout<<"\n Total Number of primary knock-on atoms: PKA = "<< N_rec<<G4endl;
  G4cout<<" Mean Kinetic Energy of PKA = "<<G4BestUnit(sum_T,"Energy")
                <<" +/- "  <<G4BestUnit(rmssum_T,"Energy")<<G4endl;
  G4cout<<" Total Frenkel-pairs Concentration: FP = " <<FP*cm3<< " #/cm3"<<G4endl;

  G4cout<<"\n----------------------------------------------------------------"<<G4endl;

  G4cout << "\n ===========================================================\n";

  //........................................................

  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

