//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: RunAction.cc,v 1.12 2004/12/02 14:53:18 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
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

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin,
                     HistoManager* histo)
:detector(det), primary(kin), ProcCounter(0), histoManager(histo)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();

  NbOfTraks0 = NbOfTraks1 = NbOfSteps0 = NbOfSteps1 = 0;
  edep = 0.;
  csdaRange = csdaRange2 = 0.;
  projRange = projRange2 = 0.;
  transvDev = transvDev2 = 0.;    
  ProcCounter = new ProcessesCount;
     
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

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  G4double dNbOfEvents = double(NbOfEvents);
    
  G4ParticleDefinition* particle = primary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();    
  G4double energy   = primary->GetParticleGun()->GetParticleEnergy();
  
  G4double length      = detector->GetSize();    
  G4Material* material = detector->GetMaterial();
  G4double density     = material->GetDensity();
   
  G4cout << "\n ======================== run summary ======================\n";
  
  G4int prec = G4cout.precision(2);
  
  G4cout << "\n The run was: " << NbOfEvents << " " << partName << " of "
         << G4BestUnit(energy,"Energy") << " through " 
	 << G4BestUnit(length,"Length") << " of "
	 << material->GetName() << " (density: " 
	 << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
	 
 G4cout << "\n ============================================================\n";
      
 G4cout.precision(3);
 
 G4cout << "\n total energy deposit: " 
        << G4BestUnit(edep/dNbOfEvents, "Energy") << G4endl;
	     
 //nb of tracks and steps per event
 //           
 G4cout << "\n nb tracks/event"
        << "   neutral: " << std::setw(10) << NbOfTraks0/dNbOfEvents
        << "   charged: " << std::setw(10) << NbOfTraks1/dNbOfEvents
        << "\n nb  steps/event"
        << "   neutral: " << std::setw(10) << NbOfSteps0/dNbOfEvents
        << "   charged: " << std::setw(10) << NbOfSteps1/dNbOfEvents
        << G4endl;
      
 //frequency of processes call       
 G4cout << "\n nb of process calls per event: \n   ";       
 for (size_t i=0; i< ProcCounter->size();i++)
     G4cout << std::setw(12) << (*ProcCounter)[i]->GetName();
           
 G4cout << "\n   ";       
 for (size_t j=0; j< ProcCounter->size();j++)
 G4cout << std::setw(12) << ((*ProcCounter)[j]->GetCounter())/dNbOfEvents;
 G4cout << G4endl;
      
 //compute csda and projected ranges, and transverse dispersion
 //
 csdaRange /= NbOfEvents; csdaRange2 /= NbOfEvents;
 G4double csdaRms = csdaRange2 - csdaRange*csdaRange;        
 if (csdaRms>0.) csdaRms = std::sqrt(csdaRms); else csdaRms = 0.;
      
 projRange /= NbOfEvents; projRange2 /= NbOfEvents;
 G4double projRms = projRange2 - projRange*projRange;        
 if (projRms>0.) projRms = std::sqrt(projRms); else projRms = 0.;
       
 transvDev /= 2*NbOfEvents; transvDev2 /= 2*NbOfEvents;
 G4double trvsRms = transvDev2 - transvDev*transvDev;        
 if (trvsRms>0.) trvsRms = std::sqrt(trvsRms); else trvsRms = 0.;
 
 //compare csda range with PhysicsTables
 //
 G4EmCalculator emCalculator;
 G4double rangeTable = 0.;
 if (particle->GetPDGCharge() != 0.)
   rangeTable = emCalculator.GetRange(energy,particle,material);
      
 G4cout << "\n---------------------------------------------------------\n";
 G4cout << " Primary particle : " ;
 G4cout << "\n CSDA Range = " << G4BestUnit(csdaRange,"Length")
        << "   rms = "        << G4BestUnit(csdaRms,  "Length");

 G4cout << "\n proj Range = " << G4BestUnit(projRange,"Length")
        << "   rms = "        << G4BestUnit(projRms,  "Length");
	     
 G4cout << "\n proj/CSDA  = " << projRange/csdaRange;
      	     
 G4cout << "\n transverse dispersion at end = " 
        << G4BestUnit(trvsRms,"Length") << G4endl;
	
 G4cout << "\n mass CSDA Range from simulation = " 
        << csdaRange*density/(g/cm2) << " g/cm2"
	<< "\n               from PhysicsTable = " 
        << rangeTable*density/(g/cm2) << " g/cm2";	
 G4cout << "\n---------------------------------------------------------\n";
                                    
  // delete and remove all contents in ProcCounter 
  while (ProcCounter->size()>0){
    OneProcessCount* aProcCount=ProcCounter->back();
    ProcCounter->pop_back();
    delete aProcCount;
  }
  delete ProcCounter;
  
  // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
 
  //save histograms      
  histoManager->save();
  
  // show Rndm status
  HepRandom::showEngineStatus(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
