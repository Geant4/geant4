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
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4ProcessTable.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronicProcess.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), fParticle(0), fEkin(0.),
  fTotalCount(0), fGammaCount(0),
  fSumTrack(0.), fSumTrack2(0.),
  fTargetXXX(false)
{
  for (G4int i=0; i<3; i++) { fPbalance[i] = 0. ; }
  for (G4int i=0; i<3; i++) { fNbGamma[i] = 0 ; }
  fPbalance[1] = DBL_MAX;
  fNbGamma[1]  = 10000;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetTargetXXX(G4bool flag)
{ 
  fTargetXXX = flag;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(G4VProcess* process) 
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

void Run::SumTrack(G4double trackl)
{
  fTotalCount++;
  fSumTrack += trackl; fSumTrack2 += trackl*trackl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountNuclearChannel(G4String name, G4double Q)
{
  std::map<G4String, NuclChannel>::iterator it = fNuclChannelMap.find(name);
  if ( it == fNuclChannelMap.end()) {
    fNuclChannelMap[name] = NuclChannel(1, Q);
  }
  else {
    NuclChannel& data = it->second;
    data.fCount++;
    data.fQ += Q;
  }       
}
                  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleCount(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMap.find(name);
  if ( it == fParticleDataMap.end()) {
    fParticleDataMap[name] = ParticleData(1, Ekin, Ekin, Ekin);
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
  }   
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Balance(G4double Pbal)
{ 
  fPbalance[0] += Pbal;
  //update min max   
  if (fTotalCount == 1) fPbalance[1] = fPbalance[2] = Pbal;  
  if (Pbal < fPbalance[1]) fPbalance[1] = Pbal;
  if (Pbal > fPbalance[2]) fPbalance[2] = Pbal;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountGamma(G4int nGamma)
{ 
  fGammaCount++;
  fNbGamma[0] += nGamma;
  //update min max   
  if (fGammaCount == 1) fNbGamma[1] = fNbGamma[2] = nGamma;  
  if (nGamma < fNbGamma[1]) fNbGamma[1] = nGamma;
  if (nGamma > fNbGamma[2]) fNbGamma[2] = nGamma;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  //primary particle info
  //
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
  
  // accumulate sums
  //
  fTotalCount   += localRun->fTotalCount;
  fGammaCount   += localRun->fGammaCount;
  fSumTrack += localRun->fSumTrack;
  fSumTrack2 += localRun->fSumTrack2;

  fPbalance[0] += localRun->fPbalance[0];
  G4double min,max;
  min = localRun->fPbalance[1]; max = localRun->fPbalance[2];
  if (fPbalance[1] > min) fPbalance[1] = min;
  if (fPbalance[2] < max) fPbalance[2] = max;

  fNbGamma[0] += localRun->fNbGamma[0];
  G4int nbmin, nbmax; 
  nbmin = localRun->fNbGamma[1]; nbmax = localRun->fNbGamma[2];
  if (fNbGamma[1] > nbmin) fNbGamma[1] = nbmin;
  if (fNbGamma[2] < nbmax) fNbGamma[2] = nbmax;
  
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
      
  //map: nuclear channels
  std::map<G4String,NuclChannel>::const_iterator itc;
  for (itc = localRun->fNuclChannelMap.begin(); 
       itc != localRun->fNuclChannelMap.end(); ++itc) {
    
    G4String name = itc->first;
    const NuclChannel& localData = itc->second;   
    if ( fNuclChannelMap.find(name) == fNuclChannelMap.end()) {
      fNuclChannelMap[name]
       = NuclChannel(localData.fCount, localData.fQ);
    }
    else {
      NuclChannel& data = fNuclChannelMap[name];   
      data.fCount += localData.fCount;
      data.fQ     += localData.fQ;
    }   
  } 
        
  //map: particles count
  std::map<G4String,ParticleData>::const_iterator itn;
  for (itn = localRun->fParticleDataMap.begin(); 
       itn != localRun->fParticleDataMap.end(); ++itn) {
    
    G4String name = itn->first;
    const ParticleData& localData = itn->second;   
    if ( fParticleDataMap.find(name) == fParticleDataMap.end()) {
      fParticleDataMap[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }
  
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun(G4bool print) 
{
  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  
  //run condition
  //
  const G4Material* material = fDetector->GetMaterial();
  G4double density = material->GetDensity();
   
  G4String Particle = fParticle->GetParticleName();    
  G4cout << "\n The run is " << numberOfEvent << " "<< Particle << " of "
         << G4BestUnit(fEkin,"Energy") << " through " 
         << G4BestUnit(fDetector->GetSize(),"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;

  if (numberOfEvent == 0) { G4cout.precision(dfprec);   return;}
             
  //frequency of processes
  //
  G4cout << "\n Process calls frequency:" << G4endl;  
  G4int survive = 0;
  std::map<G4String,G4int>::iterator it;    
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     G4cout << "\t" << procName << "= " << count;
     if (procName == "Transportation") survive = count;
  }
  G4cout << G4endl;
      
  if (survive > 0) {
    G4cout << "\n Nb of incident particles surviving after "
           << G4BestUnit(fDetector->GetSize(),"Length") << " of "
           << material->GetName() << " : " << survive << G4endl;
  }
  
  if (fTotalCount == 0) fTotalCount = 1;   //force printing anyway
  
  //compute mean free path and related quantities
  //
  G4double MeanFreePath = fSumTrack /fTotalCount;     
  G4double MeanTrack2   = fSumTrack2/fTotalCount;     
  G4double rms = std::sqrt(std::fabs(MeanTrack2 - MeanFreePath*MeanFreePath));
  G4double CrossSection = 0.0;
  if(MeanFreePath > 0.0) { CrossSection = 1./MeanFreePath; }
  G4double massicMFP = MeanFreePath*density;
  G4double massicCS  = 0.0;
  if(massicMFP > 0.0) { massicCS = 1./massicMFP; }
   
  G4cout << "\n\n MeanFreePath:\t"   << G4BestUnit(MeanFreePath,"Length")
         << " +- "                   << G4BestUnit( rms,"Length")
         << "\tmassic: "             << G4BestUnit(massicMFP, "Mass/Surface")
         << "\n CrossSection:\t"     << CrossSection*cm << " cm^-1 "
         << "\t\tmassic: "           << G4BestUnit(massicCS, "Surface/Mass")
         << G4endl;
         
  //cross section per atom (only for single material)
  //
  if (material->GetNumberOfElements() == 1) {
    G4double nbAtoms = material->GetTotNbOfAtomsPerVolume();
    G4double crossSection = CrossSection/nbAtoms;
    G4cout << " crossSection per atom:\t"
           << G4BestUnit(crossSection,"Surface") << G4endl;     
  }         
  //check cross section from G4HadronicProcessStore
  //
  G4cout << "\n Verification: "
         << "crossSections from G4HadronicProcessStore" << G4endl;
  
  G4ProcessTable* processTable  = G4ProcessTable::GetProcessTable();
  G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();
  G4double sumc1 = 0.0, sumc2 = 0.0; 
  const G4Element* element = (material->GetNumberOfElements() == 1)
    ? material->GetElement(0) : nullptr;
  for (it = fProcCounter.begin(); it != fProcCounter.end(); ++it) {
    G4String procName = it->first;
    const G4VProcess* process = processTable->FindProcess(procName, fParticle);
    PrintXS(process, material, element, store, density, sumc1, sumc2);
  }             
  if(sumc1 > 0.0) {
    G4cout << std::setw(20) << "total" << " = "
	   << G4BestUnit(sumc1, "Surface/Mass") << "\t"; 
    if(sumc2 > 0.0) { G4cout << G4BestUnit(sumc2, "Surface"); } 
    G4cout << G4endl;
  } else {
    G4cout << " not available" << G4endl;
  }
              
  //nuclear channel count
  //
  G4cout << "\n List of nuclear reactions: \n" << G4endl; 
  std::map<G4String,NuclChannel>::iterator ic;               
  for (ic = fNuclChannelMap.begin(); ic != fNuclChannelMap.end(); ic++) { 
    G4String name    = ic->first;
    NuclChannel data = ic->second;
    G4int count = data.fCount;
    G4double Q  = data.fQ/count; 
    if (print)         
      G4cout << "  " << std::setw(60) << name << ": " << std::setw(7) << count
             << "   Q = " << std::setw(wid) << G4BestUnit(Q, "Energy")
             << G4endl;           
  } 
 
  //Gamma count
  //
  if (print && (fGammaCount > 0)) {       
    G4cout << "\n" << std::setw(58) << "number of gamma or e- (ic): N = " 
           << fNbGamma[1] << " --> " << fNbGamma[2] << G4endl;
  }
 
  if (print && fTargetXXX) {
    G4cout 
      << "\n   --> NOTE: XXXX because neutronHP is unable to return target nucleus"
      << G4endl;
  }
            
  //particles count
  //
  G4cout << "\n List of generated particles:" << G4endl;
     
  std::map<G4String,ParticleData>::iterator itn;               
  for (itn = fParticleDataMap.begin(); itn != fParticleDataMap.end(); itn++) { 
    G4String name = itn->first;
    ParticleData data = itn->second;
    G4int count = data.fCount;
    G4double eMean = data.fEmean/count;
    G4double eMin = data.fEmin;
    G4double eMax = data.fEmax;    
    if (print)         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ")" << G4endl;           
  }
 
  //energy momentum balance
  //
  if (fTotalCount > 1) {
    G4double Pbmean = fPbalance[0]/fTotalCount;           
    G4cout << "\n   Momentum balance: Pmean = " 
           << std::setw(wid) << G4BestUnit(Pbmean, "Energy")
           << "\t( "  << G4BestUnit(fPbalance[1], "Energy")
           << " --> " << G4BestUnit(fPbalance[2], "Energy")
           << ") \n" << G4endl;
  }
  
  //normalize histograms      
  ////G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  ////G4double factor = 1./numberOfEvent;
  ////analysisManager->ScaleH1(3,factor);
           
  //remove all contents in fProcCounter, fCount 
  fProcCounter.clear();
  fNuclChannelMap.clear();      
  fParticleDataMap.clear();
                          
  //restore default format         
  G4cout.precision(dfprec);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::PrintXS(const G4VProcess* proc,
                  const G4Material* mat, const G4Element* elm,
		  G4HadronicProcessStore* store, G4double density,
                  G4double& sum1, G4double& sum2)
{
  if(nullptr == proc) { return; }
  G4double xs1 = store->GetCrossSectionPerVolume(fParticle, fEkin, proc, mat);
  G4double massSigma = xs1/density;
  sum1 += massSigma;
  if(nullptr != elm) {
    G4double xs2 = store->GetCrossSectionPerAtom(fParticle, fEkin, proc, elm, mat);
    sum2 += xs2;
    G4cout << "\n" << std::setw(20) << proc->GetProcessName() << " = "
	   << G4BestUnit(massSigma, "Surface/Mass") << "\t"
	   << G4BestUnit(xs2, "Surface");
  } else {
    G4cout << "\n" << std::setw(20) << proc->GetProcessName() << " = " 
	   << G4BestUnit(massSigma, "Surface/Mass");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
