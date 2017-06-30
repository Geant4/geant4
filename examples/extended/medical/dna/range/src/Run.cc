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
/// \file medical/dna/range/src/Run.cc
/// \brief Implementation of the Run class

#include "Run.hh"
#include "DetectorConstruction.hh"

#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(const DetectorConstruction* detector)
: G4Run(),
  fDetector(detector),
  fParticle(0), fEkin(0.),  
  fEdeposit(0.),  fEdeposit2(0.),
  fTrackLen(0.),  fTrackLen2(0.),
  fProjRange(0.), fProjRange2(0.),
  fPenetration(0.), fPenetration2(0.),
  fNbOfSteps(0),  fNbOfSteps2(0),
  fStepSize(0.),  fStepSize2(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary (G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin     = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEdep (G4double e)        
{
  fEdeposit  += e;
  fEdeposit2 += e*e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void Run::AddTrackLength (G4double t) 
{
  fTrackLen  += t;
  fTrackLen2 += t*t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void Run::AddProjRange (G4double x) 
{
  fProjRange  += x;
  fProjRange2 += x*x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void Run::AddPenetration (G4double x) 
{
  fPenetration  += x;
  fPenetration2 += x*x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void Run::AddStepSize (G4int nb, G4double st)
{
  fNbOfSteps  += nb; 
  fNbOfSteps2 += nb*nb;
  fStepSize   += st ; 
  fStepSize2  += st*st;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  // accumulate sums
  fEdeposit   += localRun->fEdeposit;
  fEdeposit2  += localRun->fEdeposit2;
  fTrackLen   += localRun->fTrackLen;  
  fTrackLen2  += localRun->fTrackLen2;
  fProjRange  += localRun->fProjRange; 
  fProjRange2 += localRun->fProjRange2;
  fPenetration  += localRun->fPenetration; 
  fPenetration2 += localRun->fPenetration2;
  fNbOfSteps  += localRun->fNbOfSteps ;
  fNbOfSteps2 += localRun->fNbOfSteps2;
  fStepSize   += localRun->fStepSize;  
  fStepSize2  += localRun->fStepSize2;

  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() 
{
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(2);
  
  //run conditions  
  //
  G4Material* material = fDetector->GetAbsorMaterial();
  G4double density  = material->GetDensity();       
  G4String partName = fParticle->GetParticleName();
  
  G4cout << "\n ======================= run summary ====================\n";  
  G4cout 
    << "\n The run is " << numberOfEvent << " "<< partName << " of "
    << G4BestUnit(fEkin,"Energy") << " through a sphere of radius "
    << G4BestUnit(fDetector->GetAbsorRadius(),"Length") << "of "
    << material->GetName() << " (density: " 
    << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;    

  if (numberOfEvent == 0) {
    G4cout.setf(mode,std::ios::floatfield);
    G4cout.precision(prec);  
    return;
  }
                    
  //compute track length of primary track
  //
  fTrackLen /= numberOfEvent; fTrackLen2 /= numberOfEvent;
  G4double rmsTrack = fTrackLen2 - fTrackLen*fTrackLen;      
        
  if (rmsTrack>0.) rmsTrack = std::sqrt(rmsTrack); else rmsTrack = 0.;

  G4cout.precision(3);       
  G4cout 
    << "\n Track length of primary track = " << G4BestUnit(fTrackLen,"Length")
    << " +- "                                
    << G4BestUnit( rmsTrack,"Length");
    
  //compute projected range of primary track
  //
  fProjRange /= numberOfEvent; fProjRange2 /= numberOfEvent;
  G4double rmsProj = fProjRange2 - fProjRange*fProjRange;        
  if (rmsProj>0.) rmsProj = std::sqrt(rmsProj); else rmsProj = 0.;
   
  G4cout 
    << "\n Projected range               = " 
    << G4BestUnit(fProjRange,"Length")
    << " +- "                                << G4BestUnit( rmsProj,"Length");    
    
  //compute penetration of primary track
  //
  fPenetration /= numberOfEvent; fPenetration2 /= numberOfEvent;
  G4double rmsPene = fPenetration2 - fPenetration*fPenetration;        
  if (rmsPene>0.) rmsPene = std::sqrt(rmsPene); else rmsPene = 0.;
   
  G4cout 
    << "\n Penetration                   = " 
    << G4BestUnit(fPenetration,"Length")
    << " +- "                                << G4BestUnit( rmsPene,"Length")    
    << G4endl;
    
  //

  //output file
    FILE *myFile;
    myFile = fopen ("range.txt","a");
    fprintf (myFile, "%e %e %e %e %e %e %e\n", 
     fEkin/eV, 
     fTrackLen/nm,
     rmsTrack/nm,
     fProjRange/nm,
     rmsProj/nm,
     fPenetration/nm,
     rmsPene/nm     
     );
    fclose (myFile);

  // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
    
}

