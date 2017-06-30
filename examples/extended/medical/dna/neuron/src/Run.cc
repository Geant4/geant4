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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
// $ID$
/// \file Runc.cc 
/// \brief Implementation of the Run class

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Analysis.hh"
//#include "NeuronLoadDataFile.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"  
#include "G4Molecule.hh"
#include "G4MoleculeCounter.hh"
#include "G4MoleculeGun.hh"
#include "G4H2O.hh"
#include <G4Scheduler.hh>
#include "G4MoleculeTable.hh"
#include "math.h"
#include "G4EmCalculator.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(), fDetector(det), 
  fParticle(0), fEkin(0.), 
  fLET(0.), fLET2(0.), 
  fTrackLen(0.),  fTrackLen2(0.)
{
  //fNeuronLoadParamz = new NeuronLoadDataFile(); 

  fSoma3DEdep = new G4double[fDetector->GetnbSomacomp()];
  for (G4int i=0; i<fDetector->GetnbSomacomp(); i++)
  {
   fSoma3DEdep[i]=0.;
  }  
  fDend3DEdep = new G4double[fDetector->GetnbDendritecomp()];
  for (G4int i=0; i<fDetector->GetnbDendritecomp(); i++)
  {
   fDend3DEdep[i]=0.;
  }
  fAxon3DEdep = new G4double[fDetector->GetnbAxoncomp()];
  for (G4int i=0; i<fDetector->GetnbAxoncomp(); i++)
  {
   fAxon3DEdep[i]=0.;
  }

  fEdepAll =  fEdepAll_err=fEdepMedium= fEdepMedium_err= fEdepSlice=  
  fEdepSlice_err=fEdepSoma=fEdepSoma_err=0. ;
  fEdepDend =  fEdepDend_err=fEdepAxon=  fEdepAxon_err=fEdepNeuron=  
  fEdepNeuron_err =0. ;
  fEnergyFlow   = fEnergyFlow2    = 0.;  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{
  delete[] fSoma3DEdep;
  delete[] fDend3DEdep;
  delete[] fAxon3DEdep;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddPrimaryLET(G4double let)
{ 
  fLET += let;
  fLET2 += let*let;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetTrackLength(G4double t)
{ 
  ftrackLength = t;
  fTrackLen  += t;
  fTrackLen2 += t*t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4VProcess* process) 
{
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

void Run::ParticleCount(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMap1.find(name);
  if ( it == fParticleDataMap1.end()) {
    fParticleDataMap1[name] = ParticleData(1, Ekin, Ekin, Ekin);
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

void Run::ParticleCountNeuron(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMap2.find(name);
  if ( it == fParticleDataMap2.end()) {
    fParticleDataMap2[name] = ParticleData(1, Ekin, Ekin, Ekin);
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

void Run::MoleculeCount(G4String, G4double)
{
//fMoleculeNumber = G4MoleculeCounter::Instance()
//                  ->GetNMoleculesAtTime(moleculename, Gtime);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::MoleculeCountNeuron(G4Molecule* molecule)
{
 G4String moleculename = molecule->GetName();
 std::map<G4String,G4int>::iterator it = fMoleculeNumber.find(moleculename);
  if ( it == fMoleculeNumber.end()) {
    fMoleculeNumber[moleculename] = 1;
  }
  else {
    fMoleculeNumber[moleculename]++; 
  }   
}
                       
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEflow(G4double eflow)
{ 
  fEnergyFlow += eflow;
  fEnergyFlow2 += eflow*eflow;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  //primary particle info
  //
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
  ftrackLength = localRun->ftrackLength;
  fTrackLen   += localRun->fTrackLen;  
  fTrackLen2  += localRun->fTrackLen2;  
  fLET   += localRun->fLET;  
  fLET2  += localRun->fLET2; 

  //map: processes count in simulation medium
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
  
  //map: created particles count outside neuron structure 
  std::map<G4String,ParticleData>::const_iterator itc;
  for (itc = localRun->fParticleDataMap1.begin(); 
       itc != localRun->fParticleDataMap1.end(); ++itc) {
    
    G4String name = itc->first;
    const ParticleData& localData = itc->second;   
    if ( fParticleDataMap1.find(name) == fParticleDataMap1.end()) {
      fParticleDataMap1[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap1[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }
  
  //map: created particles count inside neuron structure       
  std::map<G4String,ParticleData>::const_iterator itn;
  for (itn = localRun->fParticleDataMap2.begin(); 
       itn != localRun->fParticleDataMap2.end(); ++itn) {
    
    G4String name = itn->first;
    const ParticleData& localData = itn->second;   
    if ( fParticleDataMap2.find(name) == fParticleDataMap2.end()) {
      fParticleDataMap2[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap2[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }
  
  //map: molecule count
  //fMoleculeNumber += localRun->fMoleculeNumber;
 
  std::map<G4String,G4int>::const_iterator itm;
  for ( itm = localRun->fMoleculeNumber.begin();
        itm != localRun->fMoleculeNumber.end(); ++itm ) {

    G4String moleculeName = itm->first;
    G4int localCount = itm->second;
    if ( fMoleculeNumber.find(moleculeName) == fMoleculeNumber.end()) {
      fMoleculeNumber[moleculeName] = localCount;
    }
    else {
      fMoleculeNumber[moleculeName] += localCount;
    }  
  }  

  // hits compartments in neuron compartments
  //  
  for (G4int i=0; i<fDetector->GetnbSomacomp(); i++)
  {
   fSoma3DEdep[i] += localRun->fSoma3DEdep[i];
  }  
  for (G4int i=0; i<fDetector->GetnbDendritecomp(); i++)
  {
   fDend3DEdep[i] +=localRun->fDend3DEdep[i];
  }
  for (G4int i=0; i<fDetector->GetnbAxoncomp(); i++)
  {
   fAxon3DEdep[i] +=localRun->fAxon3DEdep[i];
  } 
  
  // accumulate sums
  //
  fEdepAll += localRun->fEdepAll;  fEdepAll_err += localRun->fEdepAll_err;
  fEdepMedium += localRun->fEdepMedium; 
  fEdepMedium_err += localRun->fEdepMedium_err;
  fEdepSlice += localRun->fEdepSlice;  
  fEdepSlice_err += localRun->fEdepSlice_err;
  fEdepSoma += localRun->fEdepSoma; fEdepSoma_err += localRun->fEdepSoma_err;
  fEdepDend += localRun->fEdepDend;  fEdepDend_err += localRun->fEdepDend_err;
  fEdepAxon += localRun->fEdepAxon;  fEdepAxon_err+= localRun->fEdepAxon_err;
  fEdepNeuron += localRun->fEdepNeuron;  
  fEdepNeuron_err += localRun->fEdepNeuron_err; 
  
  fEnergyFlow      += localRun->fEnergyFlow;
  fEnergyFlow2     += localRun->fEnergyFlow2;
  
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() 
{
  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);  
  
  // Characteristics of Primary particle 
  G4String Particle = fParticle->GetParticleName(); 
  G4double GunArea ;
  //GunArea = fParticle->GetGunArea();
  G4Material* material = fDetector->GetTargetMaterial();
  //G4double density     = material->GetDensity();   
  //compute track length of primary track
  //
  fTrackLen /= numberOfEvent; fTrackLen2 /= numberOfEvent;
  G4double rms = fTrackLen2 - fTrackLen*fTrackLen;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.; 
  
  G4int TotNbofEvents = numberOfEvent;
  G4double EdepTotal = fEdepAll;   
  G4double EdepTotal2 = fEdepAll_err;  
  EdepTotal /= TotNbofEvents; EdepTotal2 /= TotNbofEvents;
  G4double rmst = EdepTotal2 - EdepTotal*EdepTotal;
  if (rmst>0.) rmst = std::sqrt(rmst);
  else         rmst = 0.;  
  
  //Stopping Power from input Table.
  G4EmCalculator emCalculator;
  //G4double dEdxTable = 0., 
  G4double dEdxFull = 0.;
  if (fParticle->GetPDGCharge()!= 0.) { 
   // dEdxTable = emCalculator.GetDEDX(fEkin,fParticle,material);
    dEdxFull  = emCalculator.ComputeTotalDEDX(fEkin,fParticle,material);    
  }

  //Stopping Power from simulation.
  G4double meandEdx  = (EdepTotal/fTrackLen)/(keV/um); 
  G4double meandEdxerr  = (rmst/fTrackLen)/(keV/um); 
  //G4double stopPower = (meandEdx/density)/(MeV*cm2/g);    
  G4int ncols=0;
  G4int nlines = 0;
  G4float tmp, En;
  G4int Ntrav = 0;
  FILE * fp = fopen("OutputPerEvent.out","r");
  while (1)  {
    ncols = fscanf(fp," %f %f %f %f %f %f %f %f",
            &tmp, &tmp, &tmp, &tmp, &En, &tmp, &tmp, &tmp);
    if (ncols < 0) break;
    if (En>0) Ntrav++;
    nlines++;}
  fclose(fp); 
  // The surface area is calculated as spherical medium
  GunArea = fDetector->GetTotSurfMedium();  
  // Fluence dose of single paticle track
  G4double FluenceDoseperBeam = 0.160218*(dEdxFull)/(GunArea*pow(10,18)) ; 
  
 G4cout << "\n ======= The summary of simulation results 'neuron' ========\n";
 G4cout 
 << "\n  Primary particle               = " << Particle
 << "\n  Kinetic energy of beam         = " << fEkin/MeV<<" A*MeV" 
 << "\n  Particle traversals the neuron = " << Ntrav<<" of "<<numberOfEvent
 << "\n  Full LET of beam as formulas   = " <<dEdxFull/(keV/um) << " keV/um"
 << "\n  Mean LET of beam as simulation = " 
 << meandEdx << " +- "  << meandEdxerr << " keV/um"
 << "\n  Mean track length of beam      = " 
 << fTrackLen/um  << " +- " << rms << " um"  
 << "\n  Particle fluence               = " 
 << numberOfEvent/(GunArea/(cm*cm))<<" particles/cm^2"
 << "\n  Fluence dose (full)            = " 
 << numberOfEvent*FluenceDoseperBeam/(joule/kg)<<" Gy"
 << "\n  Fluence dose ber beam          = " 
 << FluenceDoseperBeam/(joule/kg) <<" Gy" << G4endl;   
 
 if (numberOfEvent == 0) { G4cout.precision(dfprec);   return;}
  
  //frequency of processes in all volume
  //
  G4cout << "\n List of generated physical process:" << G4endl;
  
  G4int index = 0;
  std::map<G4String,G4int>::iterator it;    
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     G4String space = " "; if (++index%1 == 0) space = "\n";
     G4cout << " " << std::setw(20) << procName << "="<< std::setw(7) << count
            << space;
  }
  G4cout << G4endl;

  //particles count outside neuron structure
  //
  G4cout << "\n List of generated particles outside neuron structure:" 
  << G4endl;
     
 std::map<G4String,ParticleData>::iterator itc;               
 for (itc = fParticleDataMap1.begin(); itc != fParticleDataMap1.end(); itc++) { 
    G4String name = itc->first;
    ParticleData data = itc->second;
    G4int count = data.fCount;
    G4double eMean = data.fEmean/count;
    G4double eMin = data.fEmin;
    G4double eMax = data.fEmax;    
    //-----> secondary particles flux
    G4double Eflow = data.fEmean/TotNbofEvents;  
 
    G4cout << "  " << std::setw(13) << name << " : " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ") \tEflow/event = " << G4BestUnit(Eflow, "Energy") 
     << G4endl;           
 } 
   
  //particles count inside neuron structure
 //
 G4cout << "\n Number of secondary particles inside neuron structure:" 
        << G4endl;
  
 std::map<G4String,ParticleData>::iterator itn;               
 for (itn = fParticleDataMap2.begin(); itn != fParticleDataMap2.end(); itn++) { 
    G4String name = itn->first;
    ParticleData data = itn->second;
    G4int count = data.fCount;
    //G4double eMean = data.fEmean/count;
    //G4double eMin = data.fEmin;
    //G4double eMax = data.fEmax;        
    //-----> secondary particles flux
    //G4double Eflow = data.fEmean/TotNbofEvents;        
 
    G4cout << "  " << std::setw(13) << name << " : " << std::setw(7) << count
           //<< "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           //<< "\t( "  << G4BestUnit(eMin, "Energy")
           //<< " --> " << G4BestUnit(eMax, "Energy") 
           //<< ") \tEflow/event = " << G4BestUnit(Eflow, "Energy") 
     << G4endl;
 }
 
  //molecules count inside neuron 
 // Time cut (from 1 ps to 10 ps) in class ITTrackingAction.cc 
 G4cout << "\n Number of molecular products inside neuron structure:" 
        << "\n    time: 1 ps - 10 ps "<< G4endl;
  // if (1 ns < time <= 10 ns ) MoleculeCount(molname, time) ;  
  G4int ind = 0;
  std::map<G4String,G4int>::iterator itm;    
  for (itm = fMoleculeNumber.begin(); itm != fMoleculeNumber.end(); itm++) {
     G4String moleculeName = itm->first;
     G4int    count    = itm->second;
     G4String space = " "; if (++ind%3 == 0) space = "\n";
  
     G4cout << "  " << std::setw(13) << moleculeName << " : " << std::setw(7) 
            << count << G4endl;
  }
 
  // compute total Energy and Dose deposited for all events
  //
  // EdepMedum + EdepSlice + EdepNeuron --> EdepAll
  // EdepSoma + EdepDend + EdepAxon + EdepSpines --> EdepNeuron
 G4cout << "\n Total energy (MeV) and dose (Gy) deposition :" << G4endl;
  
 G4cout << "  " 
 << std::setw(13) << "All volume:  " << std::setw(7) << fEdepAll/MeV<< " and "
 << std::setw(wid) 
 << (fEdepAll/joule)/(fDetector->GetTotMassMedium()/kg) << "\n " 
 << "  " << std::setw(13) << "Bounding slice: " 
 << std::setw(7) << (fEdepSlice+fEdepNeuron)/MeV<< " and "
 << std::setw(wid) << ((fEdepSlice+fEdepNeuron)/joule)/
    (fDetector->GetTotMassSlice()/kg) << "\n " 
 << "  " << std::setw(13) << "Neuron:   " << std::setw(7) 
 << fEdepNeuron/MeV<< " and "
 << std::setw(wid) << (fEdepNeuron/joule)/ 
    (fDetector->GetTotMassNeuron()/kg)<< "\n "  
 << "  " << std::setw(13) << "Soma:   " << std::setw(7) 
 << fEdepSoma/MeV<< " and "
 << std::setw(wid) << (fEdepSoma/joule)/ 
    (fDetector->GetMassSomaTot()/kg)<< "\n " 
 << "  " << std::setw(13) << "Dendrites:  " << std::setw(7) 
 << fEdepDend/MeV<< " and "
 << std::setw(wid) << (fEdepDend/joule)/ 
    (fDetector->GetMassDendTot()/kg)<< "\n " 
 << "  " << std::setw(13) << "Axon:   " << std::setw(7) 
 << fEdepAxon/MeV<< " and "
 << std::setw(wid) << (fEdepAxon/joule)/ 
    (fDetector->GetMassAxonTot()/kg)  
 << G4endl;
   
  // compute mean Energy and Dose deposited in hit compartments
  // 
  //G4AnalysisManager* analys = G4AnalysisManager::Instance();  
 G4int    somaCompHit = 0;
 G4double somaCompEdep = 0.;
 G4double somaCompDose = 0.;
 G4double somaCompEdep2 = 0.;
 G4double somaCompDose2 = 0.;
  // Remove old outputs  
  remove ("Soma3DEdep.out");
 for (G4int i=0; i<fDetector->GetnbSomacomp(); i++) 
 {  
  if (fSoma3DEdep[i] > 0.0) 
  {
  somaCompHit ++;
  somaCompEdep += fSoma3DEdep[i] ;
  somaCompDose += fSoma3DEdep[i]/fDetector->GetMassSomacomp(i) ;
  somaCompEdep2 += fSoma3DEdep[i]*fSoma3DEdep[i] ;
  somaCompDose2 += (fSoma3DEdep[i]/fDetector->GetMassSomacomp(i))*
     (fSoma3DEdep[i]/fDetector->GetMassSomacomp(i));
 /*G4double distSoma = 0.;
 //Fill ntuple #1
 analys->FillNtupleDColumn(1,0,i+1);
 analys->FillNtupleDColumn(1,1,distSoma);
 analys->FillNtupleDColumn(1,2,fSoma3DEdep[i]/keV);
 analys->FillNtupleDColumn(1,3,(fSoma3DEdep[i]/joule)/
         (fDetector->GetMassSomacomp(i)/kg));
 analys->AddNtupleRow(1);
 */
  
  std::ofstream WriteDataInSoma("Soma3DEdep.out", std::ios::app);
   // Index of targeted compartments 
  WriteDataInSoma //<<   i+1            << '\t' << "   " 
   // position of compartments
          <<   fDetector->GetPosSomacomp(i).x()<< '\t' << "   " 
   <<   fDetector->GetPosSomacomp(i).y()<< '\t' << "   " 
   <<   fDetector->GetPosSomacomp(i).z()<< '\t' << "   " 
   // Edep in compartments 
   <<   fSoma3DEdep[i]/keV             << '\t' << "   "  
   // Dose in compartments
   <<   (fSoma3DEdep[i]/joule)/(fDetector->GetMassSomacomp(i)/kg)
   // Dose in whole structure of Soma  
   //<<   (fSoma3DEdep[i]/joule)/(fDetector->GetMassTotSo1()/kg)
   //<< '\t' << "   " 
   << G4endl;    
  }
 }
  // compute mean Energy and Dose deposited in compartments; 
  // +- RMS : Root Mean Square  
  G4double rmsEdepS =0.;
  G4double rmsDoseS =0.;
  if (somaCompHit >0)
  {
  somaCompEdep /= somaCompHit; somaCompEdep2 /= somaCompHit;
  rmsEdepS = somaCompEdep2 - somaCompEdep*somaCompEdep;
  if (rmsEdepS>0.) rmsEdepS = std::sqrt(rmsEdepS);
  else            rmsEdepS = 0.;  
  somaCompDose /= somaCompHit; somaCompDose2 /= somaCompHit;
  rmsDoseS = somaCompDose2 - somaCompDose*somaCompDose;
  if (rmsDoseS>0.) rmsDoseS = std::sqrt(rmsDoseS);
  else            rmsDoseS = 0.;
  }
  
 G4int   DendCompHit = 0;
 G4double DendCompEdep = 0.;
 G4double DendCompDose = 0.;
 G4double DendCompEdep2 = 0.;
 G4double DendCompDose2 = 0.;
 remove ("Dend3DEdep.out");
 for (G4int i=0; i<fDetector->GetnbDendritecomp(); i++) 
 {  
  if (fDend3DEdep[i] > 0.0) 
  {
  DendCompHit ++;
  DendCompEdep += fDend3DEdep[i] ;
  DendCompDose += fDend3DEdep[i]/fDetector->GetMassDendcomp(i) ;
  DendCompEdep2 += fDend3DEdep[i]*fDend3DEdep[i] ;
  DendCompDose2 += (fDend3DEdep[i]/fDetector->GetMassDendcomp(i))*
     (fDend3DEdep[i]/fDetector->GetMassDendcomp(i)); 
 //Fill ntuple #2
 /*analys->FillNtupleDColumn(2,0,i+1);
 analys->FillNtupleDColumn(2,1,fDetector->GetDistDendsoma(i));
 analys->FillNtupleDColumn(2,2,fDend3DEdep[i]/keV);
 analys->FillNtupleDColumn(2,3,(fDend3DEdep[i]/joule)/
        (fDetector->GetMassDendcomp(i)/kg));
 analys->AddNtupleRow(2);
 */    
  
  std::ofstream WriteDataInDend("Dend3DEdep.out", std::ios::app);
  WriteDataInDend //<<   i+1            << '\t' << "   " 
   // position of compartments
   <<   fDetector->GetPosDendcomp(i).x()<< '\t' << "   " 
   <<   fDetector->GetPosDendcomp(i).y()<< '\t' << "   " 
   <<   fDetector->GetPosDendcomp(i).z()<< '\t' << "   " 
   <<   fDetector->GetDistADendSoma(i)<< '\t' << "   " 
   <<   fDetector->GetDistBDendSoma(i)<< '\t' << "   " 
   <<   fDend3DEdep[i]/keV   << '\t' << "   "  
   <<   (fDend3DEdep[i]/joule)/(fDetector->GetMassDendcomp(i)/kg)
   << G4endl;    
  }
 }
  // +- RMS : Root Mean Square  
  G4double rmsEdepD =0.;
  G4double rmsDoseD =0.;
  if (DendCompHit >0)
  {
  DendCompEdep /= DendCompHit; DendCompEdep2 /= DendCompHit;
  rmsEdepD = DendCompEdep2 - DendCompEdep*DendCompEdep;
  if (rmsEdepD>0.) rmsEdepD = std::sqrt(rmsEdepD);
  else            rmsEdepD = 0.;  
  DendCompDose /= DendCompHit; DendCompDose2 /= DendCompHit;
  rmsDoseD = DendCompDose2 - DendCompDose*DendCompDose;
  if (rmsDoseD>0.) rmsDoseD = std::sqrt(rmsDoseD);
  else            rmsDoseD = 0.;
  }
  
 G4int   AxonCompHit = 0;
 G4double AxonCompEdep = 0.;
 G4double AxonCompDose = 0.;
 G4double AxonCompEdep2 = 0.;
 G4double AxonCompDose2 = 0.;  
  remove ("Axon3DEdep.out"); 
 for (G4int i=0; i<fDetector->GetnbAxoncomp(); i++) 
 {   
  if (fAxon3DEdep[i] > 0.0) 
  {
  AxonCompHit ++;
  AxonCompEdep += fAxon3DEdep[i] ;
  AxonCompDose += fAxon3DEdep[i]/fDetector->GetMassAxoncomp(i) ;
  AxonCompEdep2 += fAxon3DEdep[i]*fAxon3DEdep[i] ;
  AxonCompDose2 += (fAxon3DEdep[i]/fDetector->GetMassAxoncomp(i))*
     (fAxon3DEdep[i]/fDetector->GetMassAxoncomp(i)); 
 //Fill ntuple #3
 /*analys->FillNtupleDColumn(3,0,i+1);
 analys->FillNtupleDColumn(3,1,fDetector->GetDistAxonsoma(i));
 analys->FillNtupleDColumn(3,2,fAxon3DEdep[i]/keV);
 analys->FillNtupleDColumn(3,3,(fAxon3DEdep[i]/joule)/
         (fDetector->GetMassAxoncomp(i)/kg));
 analys->AddNtupleRow(3);
 */     
 
  std::ofstream WriteDataInAxon("Axon3DEdep.out", std::ios::app);
  WriteDataInAxon //<<   i+1            << '\t' << "   " 
   // position of compartments
   <<   fDetector->GetPosAxoncomp(i).x()<< '\t' << "   " 
   <<   fDetector->GetPosAxoncomp(i).y()<< '\t' << "   " 
   <<   fDetector->GetPosAxoncomp(i).z()<< '\t' << "   " 
   <<   fDetector->GetDistAxonsoma(i) << '\t' << "   "  
   <<   fAxon3DEdep[i]/keV             << '\t' << "   " 
   <<   (fAxon3DEdep[i]/joule)/(fDetector->GetMassAxoncomp(i)/kg) 
   << G4endl;    
  }
 } 
  // +- RMS : Root Mean Square  
  G4double rmsEdepA =0.;
  G4double rmsDoseA =0.;
  if (AxonCompHit >0)
  {
  AxonCompEdep /= AxonCompHit; AxonCompEdep2 /= AxonCompHit;
  rmsEdepA = AxonCompEdep2 - AxonCompEdep*AxonCompEdep;
  if (rmsEdepA>0.) rmsEdepA = std::sqrt(rmsEdepA);
  else            rmsEdepA = 0.;  
  AxonCompDose /= AxonCompHit; AxonCompDose2 /= AxonCompHit;
  rmsDoseA = AxonCompDose2 - AxonCompDose*AxonCompDose;
  if (rmsDoseA>0.) rmsDoseA = std::sqrt(rmsDoseA);
  else            rmsDoseA = 0.;
  }

 G4cout << "\n Number of compartments traversed by particle tracks :" 
        << G4endl;   
 G4cout  << "  " << std::setw(13) << "Soma:  " << std::setw(7) << somaCompHit
   << " of total: "<< fDetector->GetnbSomacomp() << "\n " 
   << "  " << std::setw(13) << "Dendrites: " << std::setw(7) << DendCompHit
   << " of total: "<< fDetector->GetnbDendritecomp() <<  "\n " 
   << "  " << std::setw(13) << "Axon: " << std::setw(7) << AxonCompHit
   << " of total: "<< fDetector->GetnbAxoncomp() << "\n "    
   << G4endl;     

 G4cout << "\n Mean energy (keV) and dose (Gy) deposition in "
        << "hitting compartments :" << G4endl;   
 G4cout  
   << "  " << std::setw(13) << "Soma:  " << std::setw(7)  << somaCompEdep/keV
   << " +- "<< rmsEdepS/keV << " and "
   << somaCompDose/(joule/kg)<< " +- "<< rmsDoseS/(joule/kg)<< "\n " 
   << "  " << std::setw(13) << "Dendrites: " << std::setw(7) 
   << DendCompEdep/keV
   << " +- "<< rmsEdepD/keV << " and " 
   << DendCompDose/(joule/kg)<< " +- "<< rmsDoseD/(joule/kg)<< "\n " 
   << "  " << std::setw(13) << "Axon: " << std::setw(7) << AxonCompEdep/keV
   << " +- "<< rmsEdepA/keV << " and "
   << AxonCompDose/(joule/kg)<< " +- "<< rmsDoseA/(joule/kg)<< "\n "    
   << G4endl;  
    
  // compute mean Energy and Dose deposited per event
  // 
  fEdepNeuron /= TotNbofEvents; fEdepNeuron_err /= TotNbofEvents;
  G4double rmsEdep = fEdepNeuron_err - fEdepNeuron*fEdepNeuron;
  if (rmsEdep>0.) rmsEdep = std::sqrt(rmsEdep);
  else            rmsEdep = 0.;   
  G4double fDoseNeuron = ((fEdepNeuron/joule)/ (fDetector->
    GetTotMassNeuron()/kg)); 
  G4double fDoseNeuron_err = ((fEdepNeuron_err/joule)/ (fDetector->
    GetTotMassNeuron()/kg)); 
  G4double rmsDose = fDoseNeuron_err - fDoseNeuron*fDoseNeuron;
  if (rmsDose>0.) rmsDose = std::sqrt(rmsDose);
  else            rmsDose = 0.;     
  G4cout << "\n Mean energy (keV) and dose (Gy) deposition per event:" 
         << G4endl;  
  G4cout 
  << "  " << std::setw(13) << "Neuron:   " << std::setw(7) << fEdepNeuron/keV
  << " +- "<< rmsEdep/keV << " and "
  << std::setw(wid) << fDoseNeuron
  << " +- "<< rmsDose << "\n "   
  << G4endl;
  G4cout<< "\n Dendritic (or Axon) compartmental energy and dose deposits \n"
  << " at the distance from Soma have been written into *.out files:" 
  << "\n Dend3DEdep.out, Axon3DEdep.out, Soma3DEdep.out"
  << "\n Outputs of energy deposition per event written in data file:" 
  << "\n OutputPerEvent.out"
  << "\n " << G4endl; 
 
  //remove all contents in fProcCounter, fCount 
  fProcCounter.clear();
  fParticleDataMap2.clear();
                          
  //restore default format         
  G4cout.precision(dfprec);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
