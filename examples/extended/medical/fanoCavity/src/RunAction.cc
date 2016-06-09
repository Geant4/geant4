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
// $Id: RunAction.cc,v 1.4 2009-01-22 18:34:06 maire Exp $
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
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  CLHEP::HepRandom::showEngineStatus();
  
  //geometry
  //
  wallThickness = detector->GetWallThickness();
  wallRadius    = detector->GetWallRadius();
  mateWall      = detector->GetWallMaterial();
  densityWall   = mateWall->GetDensity();

  cavityThickness = detector->GetCavityThickness();
  cavityRadius    = detector->GetCavityRadius();
  surfaceCavity   = pi*cavityRadius*cavityRadius;
  volumeCavity    = surfaceCavity*cavityThickness;   		     
  mateCavity      = detector->GetCavityMaterial();
  densityCavity   = mateCavity->GetDensity();
  massCavity      = volumeCavity*densityCavity;
    
  //process counter
  //
  ProcCounter = new ProcessesCount;
  
  //kinetic energy of charged secondary a creation
  //
  Esecondary = Esecondary2 = 0.;
  nbSec = 0;
  
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
    
  //mean kinetic energy of secondary electrons
  //
  G4double meanEsecond = 0.;
  if (nbSec > 0) meanEsecond = Esecondary/nbSec;
  G4double rateEmean = 0.;
  // compute variation rate (%), iteration to iteration
  if (oldEmean > 0.) rateEmean = 100*(meanEsecond/oldEmean - 1.);
  oldEmean = meanEsecond;
  	        
  //beam energy fluence
  //
  G4double rBeam = wallRadius*(kinematic->GetBeamRadius());
  G4double surfaceBeam = pi*rBeam*rBeam;
  G4double beamEnergy = kinematic->GetParticleGun()->GetParticleEnergy();  
         
  //total dose in cavity
  //		   
  G4double doseCavity = EdepCavity/massCavity;
  G4double doseOverBeam = doseCavity*surfaceBeam/(NbofEvents*beamEnergy);
  G4double rateDose = 0.;
  // compute variation rate (%), iteration to iteration  
  if (oldDose > 0.) rateDose = 100*(doseOverBeam/oldDose - 1.);
  oldDose = doseOverBeam;  

  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(3);
    
  G4cout << "\n ---> NbofEvents= " << NbofEvents 
         << "   NbOfelectr= " << nbSec
	 << "   Tkin= " << G4BestUnit(meanEsecond,"Energy")
	 << " (" << rateEmean << " %)"
	 << "   NbOfelec in cav= " << PartFlowCavity[0]
	 << "   Dose/EnFluence= " << G4BestUnit(doseOverBeam,"Surface/Mass")
	 << " (" << rateDose << " %)"
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
  
  G4int NbofEvents = aRun->GetNumberOfEvent();
  if (NbofEvents == 0) return;
  
  //run conditions
  //     
  G4ParticleDefinition* particle = kinematic->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();    		         
  G4double energy = kinematic->GetParticleGun()->GetParticleEnergy();
  
  G4cout << "\n ======================== run summary ======================\n";
  
  G4int prec = G4cout.precision(3);
  
  G4cout << "\n The run consists of " << NbofEvents << " "<< partName << " of "
         << G4BestUnit(energy,"Energy") << " through 2*" 
	 << G4BestUnit(wallThickness,"Length") << " of "
	 << mateWall->GetName() << " (density: " 
	 << G4BestUnit(densityWall,"Volumic Mass") << ")" << G4endl;
	 
  G4cout << "\n the cavity is "
	 << G4BestUnit(cavityThickness,"Length") << " of "
	 << mateCavity->GetName() << " (density: " 
	 << G4BestUnit(densityCavity,"Volumic Mass") << "); Mass = " 
	 << G4BestUnit(massCavity,"Mass") << G4endl;
	 	 
  G4cout << "\n ============================================================\n";

  //frequency of processes
  //
  G4cout << "\n Process calls frequency --->";
  for (size_t i=0; i< ProcCounter->size();i++) {
     G4String procName = (*ProcCounter)[i]->GetName();
     G4int    count    = (*ProcCounter)[i]->GetCounter(); 
     G4cout << "  " << procName << "= " << count;
  }
  G4cout << G4endl;
    
  //extract cross sections with G4EmCalculator
  //
  G4EmCalculator emCalculator;  
  G4cout << "\n Gamma crossSections in wall material :";
  G4double sumc = 0.0;  
  for (size_t i=0; i< ProcCounter->size();i++) {
    G4String procName = (*ProcCounter)[i]->GetName();
    G4double massSigma = 
    emCalculator.ComputeCrossSectionPerVolume(energy,particle,
                                              procName,mateWall)/densityWall;
    if (massSigma > 0.) {
      sumc += massSigma;
      G4cout << "  " << procName << "= "
             << G4BestUnit(massSigma, "Surface/Mass");
    }	     
  }  	   
  G4cout << "   --> total= " 
         << G4BestUnit(sumc, "Surface/Mass") << G4endl;
  
  //mean kinetic energy of secondary electrons
  //
  if (nbSec == 0) return;
  G4double meanEsecond = Esecondary/nbSec, meanEsecond2 = Esecondary2/nbSec;
  G4double varianceEsec = meanEsecond2 - meanEsecond*meanEsecond;
  G4double dToverT = 0.;
  if (varianceEsec>0.) dToverT = std::sqrt(varianceEsec/nbSec)/meanEsecond;
  G4double csdaRange =
      emCalculator.GetCSDARange(meanEsecond,G4Electron::Electron(),mateWall);

  G4cout.precision(4);       
  G4cout 
    << "\n Mean energy of secondary e- = " << G4BestUnit(meanEsecond,"Energy")
    << " +- " << 100*dToverT << " %"
    << "  (--> range in wall material = "  << G4BestUnit(csdaRange,"Length")
    << ")"   
    << G4endl;
    
  //compute mass energy transfer coefficient
  //
  G4double massTransfCoef = sumc*meanEsecond/energy;
   
  G4cout << " Mass_energy_transfer coef: "  
         << G4BestUnit(massTransfCoef, "Surface/Mass")
         << G4endl;
	 
  //stopping power from EmCalculator
  //
  G4double dedxWall = 
      emCalculator.GetDEDX(meanEsecond,G4Electron::Electron(),mateWall);
  dedxWall /= densityWall;
  G4double dedxCavity = 
      emCalculator.GetDEDX(meanEsecond,G4Electron::Electron(),mateCavity);
  dedxCavity /= densityCavity;
  
  G4cout 
    << "\n StoppingPower in wall   = " 
    << G4BestUnit(dedxWall,"Energy*Surface/Mass")
    << "\n               in cavity = " 
    << G4BestUnit(dedxCavity,"Energy*Surface/Mass")
    << G4endl;  
  	
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
  	        
  //beam energy fluence
  //
  G4double rBeam = wallRadius*(kinematic->GetBeamRadius());
  G4double surfaceBeam = pi*rBeam*rBeam;
  
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
  G4double doseOverBeam = doseCavity*surfaceBeam/(NbofEvents*energy);
  
  //track length in cavity
  G4double meantrack = trkSegmCavity/PartFlowCavity[0];
  		  
  G4cout.precision(4);       
  G4cout 
    << "\n Total edep in cavity = "      << G4BestUnit(EdepCavity,"Energy")
    << " +- " << 100*dEoverE << " %"    
    << "\t Total charged trackLength = " << G4BestUnit(trkSegmCavity,"Length")
    << "   (mean value = " << G4BestUnit(meantrack,"Length") << ")"       
    << "\n Total dose in cavity = " << doseCavity/(MeV/mg) << " MeV/mg"
    << "\n Dose/EnergyFluence   = " << G4BestUnit(doseOverBeam,"Surface/Mass")
    << G4endl;
    
  //ratio simulation/theory
  //
  G4double ratio = doseOverBeam/massTransfCoef;
  G4double error = ratio*std::sqrt(dEoverE*dEoverE + dToverT*dToverT);
  
  G4cout.precision(5);  
  G4cout 
    << "\n (Dose/EnergyFluence)/Mass_energy_transfer = " << ratio 
    << " +- " << error << G4endl; 
     	 
  //compute mean step size of charged particles
  //
  stepWall /= nbStepWall; stepWall2 /= nbStepWall;
  G4double rms = stepWall2 - stepWall*stepWall;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(4);       
  G4cout 
    << "\n StepSize of ch. tracks in wall   = " 
    << G4BestUnit(stepWall,"Length") << " +- " << G4BestUnit( rms,"Length")
    << "\t (nbSteps/track = " << double(nbStepWall)/nbSec << ")";
    
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
