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
/// \file electromagnetic/TestEm15/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
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
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim)
  : G4UserRunAction(),fDetector(det), fPrimary(prim), fProcCounter(0),
    fHistoManager(0)
{
  fHistoManager = new HistoManager(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fHistoManager; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{  
  // save Rndm status
  ////G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();

  fProcCounter = new ProcessesCount;
  fTotalCount = 0;
  
  fTruePL = fTruePL2 = fGeomPL = fGeomPL2 = 0.;
  fLDispl = fLDispl2 = fPsiSpa = fPsiSpa2 = 0.;
  fTetPrj = fTetPrj2 = 0.;
  fPhiCor = fPhiCor2 = 0.;
     
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcesses(G4String procName)
{
   //does the process  already encounted ?
   size_t nbProc = fProcCounter->size();
   size_t i = 0;
   while ((i<nbProc)&&((*fProcCounter)[i]->GetName()!=procName)) i++;
   if (i == nbProc) fProcCounter->push_back( new OneProcessCount(procName));

   (*fProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  G4int  prec = G4cout.precision(5);
    
  G4Material* material = fDetector->GetMaterial();
  G4double density = material->GetDensity();
   
  G4ParticleDefinition* particle = 
                            fPrimary->GetParticleGun()->GetParticleDefinition();
  G4String Particle = particle->GetParticleName();    
  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
  G4cout << "\n The run consists of " << NbOfEvents << " "<< Particle << " of "
         << G4BestUnit(energy,"Energy") << " through " 
         << G4BestUnit(fDetector->GetBoxSize(),"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  
  //frequency of processes
  G4cout << "\n Process calls frequency --->";
  for (size_t i=0; i< fProcCounter->size();i++) {
     G4String procName = (*fProcCounter)[i]->GetName();
     G4int    count    = (*fProcCounter)[i]->GetCounter(); 
     G4cout << "\t" << procName << " = " << count;
  }
  
  if (fTotalCount > 0) {
  
    //compute path length and related quantities
    //
    G4double MeanTPL  = fTruePL /fTotalCount;     
    G4double MeanTPL2 = fTruePL2/fTotalCount;     
    G4double rmsTPL   = std::sqrt(std::fabs(MeanTPL2 - MeanTPL*MeanTPL));
    
    G4double MeanGPL  = fGeomPL /fTotalCount;     
    G4double MeanGPL2 = fGeomPL2/fTotalCount;     
    G4double rmsGPL   = std::sqrt(std::fabs(MeanGPL2 - MeanGPL*MeanGPL));
    
    G4double MeanLaD  = fLDispl /fTotalCount;     
    G4double MeanLaD2 = fLDispl2/fTotalCount;     
    G4double rmsLaD   = std::sqrt(std::fabs(MeanLaD2 - MeanLaD*MeanLaD));
    
    G4double MeanPsi  = fPsiSpa /(fTotalCount);     
    G4double MeanPsi2 = fPsiSpa2/(fTotalCount);     
    G4double rmsPsi   = std::sqrt(std::fabs(MeanPsi2 - MeanPsi*MeanPsi));
    
    G4double MeanTeta  = fTetPrj /(2*fTotalCount);     
    G4double MeanTeta2 = fTetPrj2/(2*fTotalCount);     
    G4double rmsTeta   = std::sqrt(std::fabs(MeanTeta2 - MeanTeta*MeanTeta));
    
    G4double MeanCorrel  = fPhiCor /(fTotalCount);     
    G4double MeanCorrel2 = fPhiCor2/(fTotalCount);     
    G4double rmsCorrel =
      std::sqrt(std::fabs(MeanCorrel2-MeanCorrel*MeanCorrel));
           
    G4cout << "\n\n truePathLength :\t" << G4BestUnit(MeanTPL,"Length")
           << " +- "                    << G4BestUnit( rmsTPL,"Length")
           <<   "\n geomPathLength :\t" << G4BestUnit(MeanGPL,"Length")
           << " +- "                    << G4BestUnit( rmsGPL,"Length")
           <<   "\n lateralDisplac :\t" << G4BestUnit(MeanLaD,"Length")
           << " +- "                    << G4BestUnit( rmsLaD,"Length")
           <<   "\n Psi            :\t" << MeanPsi/mrad << " mrad"
           << " +- "                    << rmsPsi /mrad << " mrad"
           <<   "  ("                   << MeanPsi/deg  << " deg"
           << " +- "                    << rmsPsi /deg  << " deg)"
           << G4endl;
    
    G4cout <<   "\n Theta_plane    :\t" << rmsTeta/mrad << " mrad"
           <<   "  ("                   << rmsTeta/deg  << " deg)"
           <<   "\n phi correlation:\t" << MeanCorrel 
           << " +- "                    << rmsCorrel
           << "  (std::cos(phi_pos - phi_dir))"                  
           << G4endl;
    
    
    //cross check from G4EmCalculator
    //
    G4cout << "\n Verification from G4EmCalculator. \n";
    
    G4EmCalculator emCal;
  
    //get transport mean free path (for multiple scattering)
    G4double MSmfp = emCal.GetMeanFreePath(energy,particle,"msc",material);
    
    //get range from restricted dedx
    G4double range = emCal.GetRangeFromRestricteDEDX(energy,particle,material);
  
    //effective facRange
    G4double efFacrange = MeanTPL/std::max(MSmfp, range);
    if (MeanTPL/range >= 0.99) efFacrange = 1.;
    
    G4cout << "\n transport mean free path :\t" << G4BestUnit(MSmfp,"Length")
           << "\n range from restrict dE/dx:\t" << G4BestUnit(range,"Length")
           << "\n ---> effective facRange  :\t" << efFacrange
           << G4endl;
    
    G4cout << "\n compute theta0 from Highland :\t"
           << ComputeMscHighland(MeanTPL)/mrad << " mrad" 
           << "  (" << ComputeMscHighland(MeanTPL)/deg << " deg)" 
           << G4endl;
                           
  } else
    G4cout<< G4endl;

  //restore default format         
  G4cout.precision(prec);         
  
  // delete and remove all contents in fProcCounter 
  while (fProcCounter->size()>0){
    OneProcessCount* aProcCount=fProcCounter->back();
    fProcCounter->pop_back();
    delete aProcCount;
  }
  delete fProcCounter;
  
  //save histograms      
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  if ( analysisManager->IsActive() ) {
  analysisManager->Write();
  analysisManager->CloseFile();
  }       

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RunAction::ComputeMscHighland(G4double pathLength)
{
 //compute the width of the Gaussian central part of the MultipleScattering
 //projected angular distribution.
 //Eur. Phys. Jour. C15 (2000) page 166, formule 23.9

 G4double t = pathLength/(fDetector->GetMaterial()->GetRadlen());
 if (t < DBL_MIN) return 0.;

 G4ParticleGun* particle = fPrimary->GetParticleGun();
 G4double T = particle->GetParticleEnergy();
 G4double M = particle->GetParticleDefinition()->GetPDGMass();
 G4double z = std::abs(particle->GetParticleDefinition()->GetPDGCharge()/eplus);

 G4double bpc = T*(T+2*M)/(T+M);
 G4double teta0 = 13.6*MeV*z*std::sqrt(t)*(1.+0.038*std::log(t))/bpc;
 return teta0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
