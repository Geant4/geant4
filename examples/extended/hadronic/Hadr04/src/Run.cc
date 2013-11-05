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
/// \file electromagnetic/TestEm11/src/Run.cc
/// \brief Implementation of the Run class
//
// $Id: Run.cc 71376 2013-06-14 07:44:50Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det),
  fNbStep1(0), fNbStep2(0),
  fTrackLen1(0.), fTrackLen2(0.),
  fTime1(0.),fTime2(0.)
{ }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleCount(G4String name, G4double Ekin)
{
  fParticleCount[name]++;
  fEmean[name] += Ekin;
  //update min max
  if (fParticleCount[name] == 1) fEmin[name] = fEmax[name] = Ekin;
  if (Ekin < fEmin[name]) fEmin[name] = Ekin;
  if (Ekin > fEmax[name]) fEmax[name] = Ekin;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumTrackLength(G4int nstep1, G4int nstep2, 
                         G4double trackl1, G4double trackl2,
                         G4double time1, G4double time2)
{
  fNbStep1   += nstep1;  fNbStep2   += nstep2;
  fTrackLen1 += trackl1; fTrackLen2 += trackl2;
  fTime1 += time1; fTime2 += time2;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ComputeStatistics() 
{
  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  
  //frequency of processes
  //
  G4cout << "\n Process calls frequency --->";  
  G4int survive = 0;
  std::map<const G4VProcess*,G4int>::iterator it;    
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first->GetProcessName();
     G4int    count    = it->second;
     G4cout << "\t" << procName << "= " << count;
     if (procName == "Transportation") survive = count;
  }
  G4cout << G4endl;
      
  if (survive > 0) {
    G4cout << "\n Nb of incident particles surviving after "
           << G4BestUnit(0.5*(fDetector->GetSize()),"Length") << " of "
           << fDetector->GetMaterial()->GetName() << " : " << survive << G4endl;
  }

 // total track length of incident neutron
 //
 G4cout << "\n Parcours of incident neutron:";
  
 G4double meanCollision1  = (G4double)fNbStep1/numberOfEvent;
 G4double meanCollision2  = (G4double)fNbStep2/numberOfEvent;
 G4double meanCollisTota  = meanCollision1 + meanCollision2;

 G4cout << "\n   nb of collisions    E>1*eV= " << meanCollision1
        << "      E<1*eV= " << meanCollision2
        << "       total= " << meanCollisTota;        
        
 G4double meanTrackLen1  = fTrackLen1/numberOfEvent;
 G4double meanTrackLen2  = fTrackLen2/numberOfEvent;
 G4double meanTrackLtot  =  meanTrackLen1 + meanTrackLen2;  

 G4cout 
   << "\n   track length        E>1*eV= " << G4BestUnit(meanTrackLen1,"Length")
   << "  E<1*eV= " << G4BestUnit(meanTrackLen2, "Length")
   << "   total= " << G4BestUnit(meanTrackLtot, "Length");   
   
 G4double meanTime1  = fTime1/numberOfEvent;
 G4double meanTime2  = fTime2/numberOfEvent;
 G4double meanTimeTo = meanTime1 + meanTime2;  

 G4cout 
   << "\n   time of flight      E>1*eV= " << G4BestUnit(meanTime1,"Time")
   << "  E<1*eV= " << G4BestUnit(meanTime2, "Time")
   << "   total= " << G4BestUnit(meanTimeTo, "Time") << G4endl;   
             
 //particles count
 //
 G4cout << "\n List of generated particles:" << G4endl;
     
 std::map<G4String,G4int>::iterator ip;               
 for (ip = fParticleCount.begin(); ip != fParticleCount.end(); ip++) { 
    G4String name = ip->first;
    G4int count   = ip->second;
    G4double eMean = fEmean[name]/count;
    G4double eMin = fEmin[name], eMax = fEmax[name];    
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ")" << G4endl;           
 }
                   
  //restore default format         
  G4cout.precision(dfprec);
           
  // remove all contents in fProcCounter 
  fProcCounter.clear();
  // remove all contents in fParticleCount
  fParticleCount.clear(); 
  fEmean.clear();  fEmin.clear(); fEmax.clear();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // accumulate sums
  //
  fNbStep1   += localRun->fNbStep1;
  fNbStep2   += localRun->fNbStep2;   
  fTrackLen1 += localRun->fTrackLen1;  
  fTrackLen2 += localRun->fTrackLen2;
  fTime1     += localRun->fTime1;  
  fTime2     += localRun->fTime2;
  
  //maps
  std::map<const G4VProcess*,G4int>::const_iterator itp;
  for ( itp = localRun->fProcCounter.begin();
        itp != localRun->fProcCounter.end(); ++itp ) {
      fProcCounter[itp->first] += itp->second;
  }    
       
  std::map<G4String,G4int>::const_iterator itn;
  for (itn = localRun->fParticleCount.begin(); 
       itn != localRun->fParticleCount.end(); ++itn) {
     fParticleCount[itn->first] += itn->second;
  }

  std::map<G4String,G4double>::const_iterator itd;
  for (itd = localRun->fEmean.begin(); 
       itd != localRun->fEmean.end(); ++itd) {
     fEmean[itd->first] += itd->second;
  }   

  for (itd = localRun->fEmin.begin(); 
       itd != localRun->fEmin.end(); ++itd) {
     G4double eminl = itd->second;
     if ( fEmin.find(itd->first) == fEmin.end() ||
          eminl < fEmin[itd->first] ) {
       fEmin[itd->first] = eminl;
     }
  }

  for (itd = localRun->fEmax.begin(); 
       itd != localRun->fEmax.end(); ++itd) {
     G4double emaxl = itd->second;
     if ( fEmax.find(itd->first) == fEmax.end() ||
          emaxl > fEmin[itd->first] ) {
       fEmax[itd->first] = emaxl;
     }
  }
  
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
