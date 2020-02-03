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
/// \file medical/fanoCavity2/src/Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Electron.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det,PrimaryGeneratorAction *kin,bool isMaster)
:G4Run(),fDetector(det), fKinematic(kin), fProcCounter(0),
  fEdepCavity(0.), fEdepCavity2(0.),
  fTrkSegmCavity(0.), fNbEventCavity(0),
  fStepWall(0.),   fStepWall2(0.),
  fStepCavity(0.), fStepCavity2(0.),
  fNbStepWall(0), fNbStepCavity(0),
  fEnergyGun(0.),  fMassWall(0.),
  fMassCavity(0.),fIsMaster(isMaster)
{

    //run conditions
    //
    G4ParticleDefinition* particleGun
                      = fKinematic->GetParticleGun()->GetParticleDefinition();
    G4String partName = particleGun->GetParticleName();
    fEnergyGun = fKinematic->GetParticleGun()->GetParticleEnergy();

    //geometry : effective wall volume
    //
    G4double cavityThickness = fDetector->GetCavityThickness();
    G4Material* mateCavity   = fDetector->GetCavityMaterial();
    G4double densityCavity   = mateCavity->GetDensity();
    fMassCavity = cavityThickness*densityCavity;

    G4double wallThickness = fDetector->GetWallThickness();
    G4Material* mateWall   = fDetector->GetWallMaterial();
    G4double densityWall   = mateWall->GetDensity();

    G4EmCalculator emCal;
    G4double RangeWall = emCal.GetCSDARange(fEnergyGun,particleGun,mateWall);
    G4double factor = 1.2;
    G4double effWallThick = factor*RangeWall;
    if ((effWallThick > wallThickness)||(effWallThick <= 0.))
      effWallThick = wallThickness;
    fMassWall = 2*effWallThick*densityWall;

    G4double massTotal     = fMassWall + fMassCavity;
    G4double fMassWallRatio = fMassWall/massTotal;
    fKinematic->RunInitialisation(effWallThick, fMassWallRatio );

    G4double massRatio = fMassCavity/fMassWall;

    //check radius
    //
    G4double worldRadius =fDetector->GetWallRadius();
    G4double RangeCavity =emCal.GetCSDARange(fEnergyGun,particleGun,mateCavity);

    //G4String partName    = particleGun->GetParticleName();


    std::ios::fmtflags mode = G4cout.flags();
    G4cout.setf(std::ios::fixed,std::ios::floatfield);
    G4int prec = G4cout.precision(3);

    G4cout << "\n ===================== run conditions =====================\n";

    G4cout << "\n The run will be " << numberOfEvent << " "<< partName << " of "
           << G4BestUnit(fEnergyGun,"Energy") << " through 2*"
           << G4BestUnit(effWallThick,"Length") << " of "
           << mateWall->GetName() << " (density: "
           << G4BestUnit(densityWall,"Volumic Mass") << "); Mass/cm2 = "
           << G4BestUnit(fMassWall*cm2,"Mass")
           << "\n csdaRange: " << G4BestUnit(RangeWall,"Length") << G4endl;

    G4cout << "\n the cavity is "
           << G4BestUnit(cavityThickness,"Length") << " of "
           << mateCavity->GetName() << " (density: "
           << G4BestUnit(densityCavity,"Volumic Mass") << "); Mass/cm2 = "
           << G4BestUnit(fMassCavity*cm2,"Mass")
           << " --> massRatio = "<< std::setprecision(6) << massRatio << G4endl;

    G4cout.precision(3);
    G4cout << " Wall radius: " << G4BestUnit(worldRadius,"Length")
           << "; range in cavity: " << G4BestUnit(RangeCavity,"Length")
           << G4endl;

    G4cout << "\n ==========================================================\n";

    //stopping power from EmCalculator
    //
    G4double dedxWall =
        emCal.GetDEDX(fEnergyGun,G4Electron::Electron(),mateWall);
    dedxWall /= densityWall;
    G4double dedxCavity =
        emCal.GetDEDX(fEnergyGun,G4Electron::Electron(),mateCavity);
    dedxCavity /= densityCavity;

    G4cout << std::setprecision(4)
           << "\n StoppingPower in wall   = "
           << G4BestUnit(dedxWall,"Energy*Surface/Mass")
           << "\n               in cavity = "
           << G4BestUnit(dedxCavity,"Energy*Surface/Mass")
           << G4endl;

    //process counter
    //
    fProcCounter = new ProcessesCount;

    //charged particles and energy flow in cavity
    //
    fPartFlowCavity[0] = fPartFlowCavity[1] = 0;
    fEnerFlowCavity[0] = fEnerFlowCavity[1] = 0.;

    //total energy deposit and charged track segment in cavity
    //
    fEdepCavity = fEdepCavity2 = fTrkSegmCavity = 0.;
    fNbEventCavity = 0;

    //stepLenth of charged particles
    //
    fStepWall = fStepWall2 = fStepCavity = fStepCavity2 =0.;
    fNbStepWall = fNbStepCavity = 0;


    // reset default formats
    G4cout.setf(mode,std::ios::floatfield);
    G4cout.precision(prec);

    //histograms
    //
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    if ( analysisManager->IsActive() ) {
      analysisManager->OpenFile();
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(G4String procName)
{
   //does the process  already encounted ?
   size_t nbProc = fProcCounter->size();
   size_t i = 0;
   while ((i<nbProc)&&((*fProcCounter)[i]->GetName()!=procName)) i++;
   if (i == nbProc) fProcCounter->push_back( new OneProcessCount(procName));

   (*fProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SurveyConvergence(G4int NbofEvents)
{
  if (NbofEvents == 0) return;




  //beam fluence
  //
  G4int Nwall   = fKinematic->GetWallCount();
  G4int Ncavity = fKinematic->GetCavityCount();
  G4double Iwall   = Nwall/fMassWall;
  G4double Icavity = Ncavity/fMassCavity;
  G4double Iratio  = Icavity/Iwall;
  G4double Itot    = NbofEvents/(fMassWall+fMassCavity);
  G4double energyFluence = fEnergyGun*Itot;

  //total dose in cavity
  //
  G4double doseCavity = fEdepCavity/fMassCavity;
  G4double ratio = doseCavity/energyFluence;
  G4double err = 100*(ratio-1.);

  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(5);

  G4cout << "--->evntNb= " << NbofEvents
         << " Nwall= " << Nwall
         << " Ncav= "  << Ncavity
         << " Ic/Iw= " << Iratio
         << " Ne-_cav= " << fPartFlowCavity[0]
         << " doseCavity/Ebeam= " << ratio
         << "  (100*(ratio-1) = " << err << " %) \n"
         << G4endl;

  // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{ // Only call by Master thread


  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(3);
  

  if (numberOfEvent == 0) return;

  //frequency of processes
  //
  G4cout << "\n Process calls frequency --->";
  for (size_t i=0; i< fProcCounter->size();i++) {
     G4String procName = (*fProcCounter)[i]->GetName();
     G4int    count    = (*fProcCounter)[i]->GetCounter(); 
     G4cout << "  " << procName << "= " << count;
  }
  G4cout << G4endl;
          
  //charged particle flow in cavity
  //
  G4cout 
    << "\n Charged particle flow in cavity :"
    << "\n      Enter --> nbParticles = " << fPartFlowCavity[0]
    << "\t Energy = " << G4BestUnit (fEnerFlowCavity[0], "Energy")
    << "\n      Exit  --> nbParticles = " << fPartFlowCavity[1]
    << "\t Energy = " << G4BestUnit (fEnerFlowCavity[1], "Energy")
    << G4endl;
             
  if (fPartFlowCavity[0] == 0) return;
                  
  G4int Nwall =  fKinematic->GetWallCount();
  G4int Ncavity = fKinematic->GetCavityCount();


  G4double Iwall   = Nwall/fMassWall;
  G4double Icavity = Ncavity/fMassCavity;
  G4double Iratio  = Icavity/Iwall;
  G4double Itot    = numberOfEvent/(fMassWall+fMassCavity);
  G4double energyFluence = fEnergyGun*Itot;  
  
  G4cout.precision(5);       
  G4cout 
    << "\n beamFluence in wall = " << Nwall
    << "\t in cavity = " << Ncavity
    << "\t Icav/Iwall = " << Iratio        
    << "\t energyFluence = " << energyFluence/(MeV*cm2/mg) << " MeV*cm2/mg"
    << G4endl;
  
  //error on Edep in cavity
  //
  if (fNbEventCavity == 0) return;
  G4double meanEdep  = fEdepCavity/fNbEventCavity;
  G4double meanEdep2 = fEdepCavity2/fNbEventCavity;
  G4double varianceEdep = meanEdep2 - meanEdep*meanEdep;
  G4double dEoverE = 0.;
  if(varianceEdep>0.) dEoverE = std::sqrt(varianceEdep/fNbEventCavity)/meanEdep;
               
  //total dose in cavity
  //                   
  G4double doseCavity = fEdepCavity/fMassCavity;

  G4double ratio = doseCavity/energyFluence, error = ratio*dEoverE;
                    
  G4cout 
    << "\n Total edep in cavity = "      << G4BestUnit(fEdepCavity,"Energy")
    << " +- " << 100*dEoverE << " %"        
    << "\n Total dose in cavity = " << doseCavity/(MeV*cm2/mg) << " MeV*cm2/mg"
    << " +- " << 100*dEoverE << " %"          
    << "\n\n DoseCavity/EnergyFluence = " << ratio 
    << " +- " << error << G4endl;
    

  //track length in cavity
  G4double meantrack = fTrkSegmCavity/fPartFlowCavity[0];
  
  G4cout.precision(4); 
  G4cout  
    << "\n Total charged trackLength in cavity = " 
    << G4BestUnit(fTrkSegmCavity,"Length")
    << "   (mean value = " << G4BestUnit(meantrack,"Length") << ")"       
    << G4endl;
                  
  //compute mean step size of charged particles
  //
  fStepWall /= fNbStepWall; fStepWall2 /= fNbStepWall;
  G4double rms = fStepWall2 - fStepWall*fStepWall;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;
  G4double nbTrackWall = fKinematic->GetWallCount();

  G4cout 
    << "\n StepSize of ch. tracks in wall   = " 
    << G4BestUnit(fStepWall,"Length") << " +- " << G4BestUnit( rms,"Length")
    << "\t (nbSteps/track = " << double(fNbStepWall)/nbTrackWall << ")";
    
  fStepCavity /= fNbStepCavity; fStepCavity2 /= fNbStepCavity;
  rms = fStepCavity2 - fStepCavity*fStepCavity;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout 
    << "\n StepSize of ch. tracks in cavity = " 
    << G4BestUnit(fStepCavity,"Length") << " +- " << G4BestUnit( rms,"Length")
    << "\t (nbSteps/track = " <<double(fNbStepCavity)/fPartFlowCavity[0] << ")";
        
  G4cout << G4endl;
  
   // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
  
  // delete and remove all contents in fProcCounter 
  while (fProcCounter->size()>0){
    OneProcessCount* aProcCount=fProcCounter->back();
    fProcCounter->pop_back();
    delete aProcCount;
  }
  delete fProcCounter;
  
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run) {

   const Run* localRun = static_cast<const Run*>(run);

    // Merge Run variables
      fPartFlowCavity[0]+= localRun->fPartFlowCavity[0];
      fPartFlowCavity[1]+= localRun->fPartFlowCavity[1];
      fEnerFlowCavity[0]+= localRun->fEnerFlowCavity[0];
      fEnerFlowCavity[1]+= localRun->fEnerFlowCavity[1];
      fEdepCavity     += localRun->fEdepCavity;
      fEdepCavity2      += localRun->fEdepCavity2;
      fTrkSegmCavity    += localRun->fTrkSegmCavity;
      fNbEventCavity    += localRun->fNbEventCavity;
      fStepWall       += localRun->fStepWall;
      fStepWall2      += localRun->fStepWall2;
    fStepCavity     += localRun->fStepCavity;
    fStepCavity2    += localRun->fStepCavity2;
    fNbStepWall       += localRun->fNbStepWall;
    fNbStepCavity     += localRun->fNbStepCavity;

    // Merge PrimaryGenerator variables
    fKinematic->AddWallCount(localRun->fKinematic->GetWallCount());
    fKinematic->AddCavityCount(localRun->fKinematic->GetCavityCount());

    // Merge ProcessCount varaibles
    std::vector<OneProcessCount*>::iterator it;
    for (it = localRun->fProcCounter->begin();
         it !=localRun->fProcCounter->end(); it++  )
    {
      OneProcessCount* process = *it;
      for ( G4int i = 0; i < process->GetCounter() ; i++)
        this->CountProcesses(process->GetName());
    }

  G4Run::Merge(run);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
