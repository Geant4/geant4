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
// $Id: RunAction.cc,v 1.3 2007-11-05 13:19:16 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Electron.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin, 
                     HistoManager* histo)
:detector(det),kinematic(kin),ProcCounter(0),histoManager(histo)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{    
  // do not save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  CLHEP::HepRandom::showEngineStatus();
  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
    
  G4int NbofEvents = aRun->GetNumberOfEventToBeProcessed();
  if (NbofEvents == 0) return;
  
  //run conditions
  //     
  G4ParticleDefinition* particleGun 
                    = kinematic->GetParticleGun()->GetParticleDefinition();
  G4String partName = particleGun->GetParticleName();    		         
  energyGun = kinematic->GetParticleGun()->GetParticleEnergy();
  
  //geometry : effective wall volume
  //  
  G4double cavityThickness = detector->GetCavityThickness();	     
  G4Material* mateCavity   = detector->GetCavityMaterial();
  G4double densityCavity   = mateCavity->GetDensity();
  massCavity = cavityThickness*densityCavity;
      
  G4double wallThickness = detector->GetWallThickness();
  G4Material* mateWall   = detector->GetWallMaterial();
  G4double densityWall   = mateWall->GetDensity();
  
  G4EmCalculator emCal;
  G4double RangeWall = emCal.GetCSDARange(energyGun,particleGun,mateWall);
  G4double factor = 1.2;
  G4double effWallThick = factor*RangeWall;
  if ((effWallThick > wallThickness)||(effWallThick <= 0.))
    effWallThick = wallThickness;
  massWall = 2*effWallThick*densityWall;  
  
  G4double massTotal     = massWall + massCavity;
  G4double massWallRatio = massWall/massTotal;  
  kinematic->RunInitialisation(effWallThick, massWallRatio ); 
     
  G4double massRatio = massCavity/massWall;
  
  //check radius
  //
  G4double worldRadius = detector->GetWorldRadius();
  G4double RangeCavity = emCal.GetCSDARange(energyGun,particleGun,mateCavity);
    
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(3);
  
  G4cout << "\n ======================== run conditions =====================\n";
  
  G4cout << "\n The run will be " << NbofEvents << " "<< partName << " of "
         << G4BestUnit(energyGun,"Energy") << " through 2*" 
	 << G4BestUnit(effWallThick,"Length") << " of "
	 << mateWall->GetName() << " (density: " 
	 << G4BestUnit(densityWall,"Volumic Mass") << "); Mass/cm2 = "
	 << G4BestUnit(massWall*cm2,"Mass") 
	 << "\n csdaRange: " << G4BestUnit(RangeWall,"Length") << G4endl;
	 
  G4cout << "\n the cavity is "
	 << G4BestUnit(cavityThickness,"Length") << " of "
	 << mateCavity->GetName() << " (density: " 
	 << G4BestUnit(densityCavity,"Volumic Mass") << "); Mass/cm2 = " 
	 << G4BestUnit(massCavity*cm2,"Mass") 
	 << " --> massRatio = " << std::setprecision(6) << massRatio << G4endl;
	  
  G4cout.precision(3);	 
  G4cout << " World radius: " << G4BestUnit(worldRadius,"Length")
         << "; range in cavity: " << G4BestUnit(RangeCavity,"Length")
	 << G4endl;	 
	 	 
  G4cout << "\n ============================================================\n";
      	 
  //stopping power from EmCalculator
  //
  G4double dedxWall = 
      emCal.GetDEDX(energyGun,G4Electron::Electron(),mateWall);
  dedxWall /= densityWall;
  G4double dedxCavity = 
      emCal.GetDEDX(energyGun,G4Electron::Electron(),mateCavity);
  dedxCavity /= densityCavity;
  
  G4cout << std::setprecision(4)
         << "\n StoppingPower in wall   = " 
         << G4BestUnit(dedxWall,"Energy*Surface/Mass")
         << "\n               in cavity = " 
         << G4BestUnit(dedxCavity,"Energy*Surface/Mass")
         << G4endl;
                  
  //process counter
  //
  ProcCounter = new ProcessesCount;
  
  //charged particles and energy flow in cavity
  //
  PartFlowCavity[0] = PartFlowCavity[1] = 0;
  EnerFlowCavity[0] = EnerFlowCavity[1] = 0.;
       
  //total energy deposit and charged track segment in cavity
  //
  EdepCavity = EdepCavity2 = trkSegmCavity = 0.;
  nbEventCavity = 0;
      
  //stepLenth of charged particles
  //
  stepWall = stepWall2 = stepCavity = stepCavity2 =0.;
  nbStepWall = nbStepCavity = 0;
    
  //histograms
  //
  histoManager->book();
  
  // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcesses(G4String procName)
{
   //does the process  already encounted ?
   size_t nbProc = ProcCounter->size();
   size_t i = 0;
   while ((i<nbProc)&&((*ProcCounter)[i]->GetName()!=procName)) i++;
   if (i == nbProc) ProcCounter->push_back( new OneProcessCount(procName));

   (*ProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SurveyConvergence(G4int NbofEvents)
{  
  if (NbofEvents == 0) return;
    
  	        
  //beam fluence
  //
  G4int Nwall   = kinematic->GetWallCount();
  G4int Ncavity = kinematic->GetCavityCount();
  G4double Iwall   = Nwall/massWall;    
  G4double Icavity = Ncavity/massCavity;
  G4double Iratio  = Icavity/Iwall;
  G4double Itot    = NbofEvents/(massWall+massCavity);
  G4double energyFluence = energyGun*Itot;
           
  //total dose in cavity
  //		   
  G4double doseCavity = EdepCavity/massCavity;
  G4double ratio = doseCavity/energyFluence;
  G4double err = 100*(ratio-1.);

  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(5);
  
  G4cout << "\n--->evntNb= " << NbofEvents 
         << " Nwall= " << Nwall
         << " Ncav= "  << Ncavity
         << " Ic/Iw= " << Iratio        
	 << " Ne-_cav= " << PartFlowCavity[0]
	 << " doseCavity/Ebeam= " << ratio 
	 << "  (100*(ratio-1) = " << err << " %)"
	 << G4endl;
	 
  // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(3);
  
  G4int NbofEvents = aRun->GetNumberOfEvent();
  if (NbofEvents == 0) return;

  //frequency of processes
  //
  G4cout << "\n Process calls frequency --->";
  for (size_t i=0; i< ProcCounter->size();i++) {
     G4String procName = (*ProcCounter)[i]->GetName();
     G4int    count    = (*ProcCounter)[i]->GetCounter(); 
     G4cout << "  " << procName << "= " << count;
  }
  G4cout << G4endl;
  	
  //charged particle flow in cavity
  //
  G4cout 
    << "\n Charged particle flow in cavity :"
    << "\n      Enter --> nbParticles = " << PartFlowCavity[0]
    << "\t Energy = " << G4BestUnit (EnerFlowCavity[0], "Energy")
    << "\n      Exit  --> nbParticles = " << PartFlowCavity[1]
    << "\t Energy = " << G4BestUnit (EnerFlowCavity[1], "Energy")
    << G4endl;
             
  if (PartFlowCavity[0] == 0) return;
  	        
  //beam fluence
  //
  G4int Nwall   = kinematic->GetWallCount();
  G4int Ncavity = kinematic->GetCavityCount();  
  G4double Iwall   = Nwall/massWall;
  G4double Icavity = Ncavity/massCavity;
  G4double Iratio  = Icavity/Iwall;
  G4double Itot    = NbofEvents/(massWall+massCavity);
  G4double energyFluence = energyGun*Itot;  
  
  G4cout.precision(5);       
  G4cout 
    << "\n beamFluence in wall = " << Nwall
    << "\t in cavity = " << Ncavity
    << "\t Icav/Iwall = " << Iratio        
    << "\t energyFluence = " << energyFluence/(MeV*cm2/mg) << " MeV*cm2/mg"
    << G4endl;
  
  //error on Edep in cavity
  //
  if (nbEventCavity == 0) return;
  G4double meanEdep  = EdepCavity/nbEventCavity;
  G4double meanEdep2 = EdepCavity2/nbEventCavity;
  G4double varianceEdep = meanEdep2 - meanEdep*meanEdep;
  G4double dEoverE = 0.;
  if(varianceEdep>0.) dEoverE = std::sqrt(varianceEdep/nbEventCavity)/meanEdep;
               
  //total dose in cavity
  //		   
  G4double doseCavity = EdepCavity/massCavity;
  G4double ratio = doseCavity/energyFluence, error = ratio*dEoverE;
  		  
  G4cout 
    << "\n Total edep in cavity = "      << G4BestUnit(EdepCavity,"Energy")
    << " +- " << 100*dEoverE << " %"        
    << "\n Total dose in cavity = " << doseCavity/(MeV*cm2/mg) << " MeV*cm2/mg"
    << " +- " << 100*dEoverE << " %"          
    << "\n\n DoseCavity/EnergyFluence = " << ratio 
    << " +- " << error << G4endl;
    

  //track length in cavity
  G4double meantrack = trkSegmCavity/PartFlowCavity[0];
  
  G4cout.precision(4); 
  G4cout  
    << "\n Total charged trackLength in cavity = " 
    << G4BestUnit(trkSegmCavity,"Length")
    << "   (mean value = " << G4BestUnit(meantrack,"Length") << ")"       
    << G4endl;
         	 
  //compute mean step size of charged particles
  //
  stepWall /= nbStepWall; stepWall2 /= nbStepWall;
  G4double rms = stepWall2 - stepWall*stepWall;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;
  G4double nbTrackWall = kinematic->GetWallCount();

  G4cout 
    << "\n StepSize of ch. tracks in wall   = " 
    << G4BestUnit(stepWall,"Length") << " +- " << G4BestUnit( rms,"Length")
    << "\t (nbSteps/track = " << double(nbStepWall)/nbTrackWall << ")";
    
  stepCavity /= nbStepCavity; stepCavity2 /= nbStepCavity;
  rms = stepCavity2 - stepCavity*stepCavity;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout 
    << "\n StepSize of ch. tracks in cavity = " 
    << G4BestUnit(stepCavity,"Length") << " +- " << G4BestUnit( rms,"Length")
    << "\t (nbSteps/track = " << double(nbStepCavity)/PartFlowCavity[0] << ")";
        
  G4cout << G4endl;
  
   // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
  
  // delete and remove all contents in ProcCounter 
  while (ProcCounter->size()>0){
    OneProcessCount* aProcCount=ProcCounter->back();
    ProcCounter->pop_back();
    delete aProcCount;
  }
  delete ProcCounter;
  
  // save histograms
  histoManager->save();
 
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
