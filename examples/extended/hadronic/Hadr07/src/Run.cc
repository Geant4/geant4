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
/// \file Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"

#include "EventAction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Material.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* detector)
: G4Run(),
  fDetector(detector),
  fParticle(0), fEkin(0.)
{
  for (G4int i=0; i<3; ++i) { fStatus[i] = 0; fTotEdep[i] = 0.; }
  fTotEdep[1] = joule;
  for (G4int i=0; i<kMaxAbsor; ++i) {
    fEdeposit[i] = 0.; fEmin[i] = joule; fEmax[i] = 0.;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary (G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin     = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4VProcess* process) 
{
  if (process == nullptr) return;
  G4String procName = process->GetProcessName();
  std::map<G4String,G4int>::iterator it = fProcCounter.find(procName);
  if ( it == fProcCounter.end()) {
    fProcCounter[procName] = 1;
  }
  else {
    fProcCounter[procName]++; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleCount(G4int k, G4String name, G4double Ekin, G4double meanLife)
{
 std::map<G4String, ParticleData>::iterator it = fParticleDataMap[k].find(name);
  if ( it == fParticleDataMap[k].end()) {
    (fParticleDataMap[k])[name] = ParticleData(1, Ekin, Ekin, Ekin, meanLife);
  }
  else {
    ParticleData& data = it->second;
    data.fCount++;
    data.fEmean += Ekin;
    //update min max
    G4double emin = data.fEmin;
    if (Ekin < emin) data.fEmin = Ekin;
    G4double emax = data.fEmax;
    if (Ekin > emax) data.fEmax = Ekin;
    data.fTmean = meanLife; 
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEdep (G4int i, G4double e)        
{
  if (e > 0.) {
    fEdeposit[i]  += e;
    if (e < fEmin[i]) fEmin[i] = e;
    if (e > fEmax[i]) fEmax[i] = e;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddTotEdep (G4double e)        
{
  if (e > 0.) {
    fTotEdep[0]  += e;
    if (e < fTotEdep[1]) fTotEdep[1] = e;
    if (e > fTotEdep[2]) fTotEdep[2] = e;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
      
void Run::AddTrackStatus (G4int i)    
{
  fStatus[i]++ ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  // Edep in absorbers
  //
  G4int nbOfAbsor = fDetector->GetNbOfAbsor();
  for (G4int i=1; i<=nbOfAbsor; ++i) {
    fEdeposit[i]  += localRun->fEdeposit[i];
    // min, max
    G4double min,max;
    min = localRun->fEmin[i]; max = localRun->fEmax[i];
    if (fEmin[i] > min) fEmin[i] = min;
    if (fEmax[i] < max) fEmax[i] = max;
  }
   
  for (G4int i=0; i<3; ++i)  fStatus[i] += localRun->fStatus[i];
  
  // total Edep
  fTotEdep[0] += localRun->fTotEdep[0];
  G4double min,max;
  min = localRun->fTotEdep[1]; max = localRun->fTotEdep[2];
  if (fTotEdep[1] > min) fTotEdep[1] = min;
  if (fTotEdep[2] < max) fTotEdep[2] = max;

  //map: processes count
  std::map<G4String,G4int>::const_iterator itp;
  for ( itp = localRun->fProcCounter.begin();
        itp != localRun->fProcCounter.end(); ++itp ) {

    G4String procName = itp->first;
    G4int localCount = itp->second;
    if ( fProcCounter.find(procName) == fProcCounter.end()) {
      fProcCounter[procName] = localCount;
    }
    else {
      fProcCounter[procName] += localCount;
    }  
  }
  
  //map: created particles in absorbers count
  for (G4int k=0; k<=nbOfAbsor; ++k) {
    std::map<G4String,ParticleData>::const_iterator itc;
    for (itc = localRun->fParticleDataMap[k].begin(); 
         itc != localRun->fParticleDataMap[k].end(); ++itc) {

      G4String name = itc->first;
      const ParticleData& localData = itc->second;   
      if ( fParticleDataMap[k].find(name) == fParticleDataMap[k].end()) {
        (fParticleDataMap[k])[name]
         = ParticleData(localData.fCount, 
                        localData.fEmean, 
                        localData.fEmin, 
                        localData.fEmax,
			localData.fTmean);
      }
      else {
        ParticleData& data = (fParticleDataMap[k])[name];   
        data.fCount += localData.fCount;
        data.fEmean += localData.fEmean;
        G4double emin = localData.fEmin;
        if (emin < data.fEmin) data.fEmin = emin;
        G4double emax = localData.fEmax;
        if (emax > data.fEmax) data.fEmax = emax;
	data.fTmean = localData.fTmean; 
      }
    }
  }

  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() 
{
    G4int prec = 5, wid = prec + 2;  
    G4int dfprec = G4cout.precision(prec);

  //run conditions
  //     
  G4String partName = fParticle->GetParticleName();
  G4int nbOfAbsor   = fDetector->GetNbOfAbsor();
  
  G4cout << "\n ======================== run summary =====================\n";
  G4cout 
    << "\n The run is " << numberOfEvent << " "<< partName << " of "
    << G4BestUnit(fEkin,"Energy") 
    << " through "  << nbOfAbsor << " absorbers: \n";
  for (G4int i=1; i<= nbOfAbsor; i++) {
     G4Material* material = fDetector->GetAbsorMaterial(i);
     G4double thickness = fDetector->GetAbsorThickness(i);
     G4double density = material->GetDensity();
     G4cout << std::setw(5) << i
            << std::setw(10) << G4BestUnit(thickness,"Length") << " of "
            << material->GetName() << " (density: " 
            << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  }         

  if (numberOfEvent == 0) { G4cout.precision(dfprec);  return;}
  
  G4cout.precision(3);
  
  //frequency of processes
  //
  G4cout << "\n Process calls frequency :" << G4endl;
  G4int index = 0;
  std::map<G4String,G4int>::iterator it;    
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     G4String space = " "; if (++index%3 == 0) space = "\n";
     G4cout << " " << std::setw(20) << procName << "="<< std::setw(7) << count
            << space;
  }
  G4cout << G4endl;
  
  //Edep in absorbers
  //
  for (G4int i=1; i<= nbOfAbsor; i++) {
     fEdeposit[i] /= numberOfEvent;

     G4cout 
       << "\n Edep in absorber " << i << " = " 
       << G4BestUnit(fEdeposit[i],"Energy")
       << "\t(" << G4BestUnit(fEmin[i], "Energy")
       << "-->" << G4BestUnit(fEmax[i], "Energy")
       << ")";
  }
  G4cout << G4endl;

  if (nbOfAbsor > 1) {
    fTotEdep[0] /= numberOfEvent;
    G4cout 
      << "\n Edep in all absorbers = " << G4BestUnit(fTotEdep[0],"Energy")
      << "\t(" << G4BestUnit(fTotEdep[1], "Energy")
      << "-->" << G4BestUnit(fTotEdep[2], "Energy")
      << ")" << G4endl;
  }

  //particles count in absorbers
  //
  for (G4int k=1; k<= nbOfAbsor; k++) {
  G4cout << "\n List of created particles in absorber " 
         << k << ":" << G4endl;

    std::map<G4String,ParticleData>::iterator itc;               
    for (itc  = fParticleDataMap[k].begin();
         itc != fParticleDataMap[k].end(); itc++) {
       G4String name = itc->first;
       ParticleData data = itc->second;
       G4int count = data.fCount;
       G4double eMean = data.fEmean/count;
       G4double eMin = data.fEmin;
       G4double eMax = data.fEmax;
       G4double meanLife = data.fTmean;    

       G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
              << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
              << "\t( "  << G4BestUnit(eMin, "Energy")
              << " --> " << G4BestUnit(eMax, "Energy") << ")";
       if (meanLife >= 0.)
         G4cout << "\tmean life = " << G4BestUnit(meanLife, "Time") << G4endl;
       else G4cout << "\tstable" << G4endl;	            
    }
  }
  //particles emerging from absorbers
  //
  G4cout << "\n List of particles emerging from absorbers :" << G4endl;
  
  std::map<G4String,ParticleData>::iterator itc;
  for (itc  = fParticleDataMap[0].begin();
       itc != fParticleDataMap[0].end(); itc++) {
    G4String name = itc->first;
    ParticleData data = itc->second;
    G4int count = data.fCount;
    G4double eMean = data.fEmean/count;
    G4double eMin = data.fEmin;
    G4double eMax = data.fEmax;
    ///G4double meanLife = data.fTmean;    

    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ")" << G4endl;
  }

  //transmission coefficients
  //
  G4double dNofEvents = double(numberOfEvent);
  G4double absorbed  = 100.*fStatus[0]/dNofEvents;
  G4double transmit  = 100.*fStatus[1]/dNofEvents;
  G4double reflected = 100.*fStatus[2]/dNofEvents;  

  G4cout.precision(2);       
  G4cout 
    << "\n Nb of events with primary absorbed = "  << absorbed  << " %,"
    << "   transmit = "  << transmit  << " %,"
    << "   reflected = " << reflected << " %" << G4endl;

  // normalize histograms of longitudinal energy profile
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int ih = 10;
  G4double binWidth = analysisManager->GetH1Width(ih)
                     *analysisManager->GetH1Unit(ih);
  G4double fac = (1./(numberOfEvent*binWidth))*(mm/MeV);
  analysisManager->ScaleH1(ih,fac);

  //remove all contents in fProcCounter, fCount 
  fProcCounter.clear();
  for (G4int k=0; k<= nbOfAbsor; k++) fParticleDataMap[k].clear();

  // reset default formats
  G4cout.precision(dfprec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
