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
/// \file medical/GammaTherapy/src/Histo.cc
/// \brief Implementation of the Histo class
//
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   Histo
//
//
// Author:      V.Ivanchenko 30/01/01
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Histo.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4RootAnalysisManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo* Histo::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo* Histo::GetPointer()
{
  if(!fManager) {
    static Histo manager;
    fManager = &manager;
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::Histo()
{
  fVerbose   = 1;
  fHistName  = "";
  fNHisto    = 10;
  fMaxEnergy = 50.0*MeV;
  fNBinsZ    = 60;
  fNBinsR    = 80;
  fNBinsE    = 200;
  fScoreBin  = 0;

  fAbsorberZ = 300.*mm;
  fAbsorberR = 200.*mm;
  fScoreZ    = 100.*mm;

  fGamma    = G4Gamma::Gamma();
  fElectron = G4Electron::Electron();
  fPositron = G4Positron::Positron();

  fCheckVolume = fGasVolume = fPhantom = fTarget1 = fTarget2 = 0;

  fAnalysisManager = 0;

  fStepZ = fStepR = fStepE = fNormZ = fSumR = 0.0;
  fNevt = fNelec = fNposit = fNgam = fNstep = fNgamPh = 
    fNgamTar = fNeTar = fNePh = fNstepTarget = 0;

  fHisto.resize(fNHisto, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::~Histo()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::BeginOfHisto()
{
  G4cout << "### Histo start initialisation nHisto= " << fNHisto << G4endl;

  fNevt = fNelec = fNposit= fNgam = fNstep = fNgamPh = fNgamTar = 
    fNeTar = fNePh = fNstepTarget = 0;
  fSumR = 0.0;
  if(fNBinsR>1000) { SetNumberDivR(40); }

  fStepZ = fAbsorberZ/(G4double)fNBinsZ;
  fStepR = fAbsorberR/(G4double)fNBinsR;
  fStepE = fMaxEnergy/(G4double)fNBinsE;
  fScoreBin = (G4int)(fScoreZ/fStepZ + 0.5);

  G4cout << "   "<< fNBinsR << " bins R   stepR= " << fStepR/mm << " mm " 
         << G4endl;
  G4cout << "   "<< fNBinsZ << " bins Z   stepZ= " << fStepZ/mm << " mm " 
         << G4endl;
  G4cout << "   "<< fNBinsE << " bins E   stepE= " << fStepE/MeV << " MeV " 
         << G4endl;
  G4cout << "   "<< fScoreBin << "-th bin in Z is used for R distribution" 
         << G4endl;

  fVolumeR.resize(fNBinsR);
  fEdep.resize(fNBinsR, 0.0);
  fGammaE.resize(fNBinsE, 0.0);

  G4double r1 = 0.0;
  G4double r2 = fStepR;
  for(G4int i=0; i<fNBinsR; ++i) {
    fVolumeR[i] = cm*cm/(pi*(r2*r2 - r1*r1));
    r1 = r2;
    r2 += fStepR;
  }

  if(fHistName != "") {
    if(!fAnalysisManager) {
      fAnalysisManager = G4RootAnalysisManager::Instance(); 
      BookHisto();
    }
    if(fVerbose > 0) {
      G4cout << "Histo: Histograms are booked and run has been started"
             << G4endl;
    }
  } else if(fVerbose > 0) {
    G4cout << "Histo: Histograms are not booked because file name is not set"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::EndOfHisto()
{

  G4cout << "Histo: End of run actions are started" << G4endl;

  // average

  G4cout<<"========================================================"<<G4endl;
  G4double x = (G4double)fNevt;
  if(fNevt > 0) { x = 1.0/x; }
  G4double xe = x*(G4double)fNelec;
  G4double xg = x*(G4double)fNgam;
  G4double xp = x*(G4double)fNposit;
  G4double xs = x*(G4double)fNstep;
  G4double xph= x*(G4double)fNgamPh;
  G4double xes= x*(G4double)fNstepTarget;
  G4double xgt= x*(G4double)fNgamTar;
  G4double xet= x*(G4double)fNeTar;
  G4double xphe= x*(G4double)fNePh;

  G4cout                    << "Number of events                             " 
                            << std::setprecision(8) << fNevt <<G4endl;
  G4cout << std::setprecision(4) << "Average number of e-                         " 
         << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma                      " 
         << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+                         " 
         << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps in the phantom       " 
         << xs << G4endl;
  G4cout << std::setprecision(4) << "Average number of e- steps in the target     " 
         << xes << G4endl;
  G4cout << std::setprecision(4) << "Average number of g  produced in the target  " 
         << xgt << G4endl;
  G4cout << std::setprecision(4) << "Average number of e- produced in the target  " 
         << xet << G4endl;
  G4cout << std::setprecision(4) << "Average number of g produced in the phantom  " 
         << xph << G4endl;
  G4cout << std::setprecision(4) << "Average number of e- produced in the phantom " 
         << xphe << G4endl;
  G4cout << std::setprecision(4) << "Total fGamma fluence in front of the phantom " 
         << x*fSumR/MeV
         << " MeV " << G4endl;
  G4cout<<"========================================================"<<G4endl;
  G4cout<<G4endl;
  G4cout<<G4endl;

  G4double sum = 0.0;
  for(G4int i=0; i<fNBinsR; i++) { 
    fEdep[i] *= x; 
    sum += fEdep[i];
  }

  if(fAnalysisManager) {

    // normalise histograms
    for(G4int i=0; i<fNHisto; i++) {
      fAnalysisManager->GetH1(fHisto[i])->scale(x);
    }
    G4double nr = fEdep[0];
    if(nr > 0.0) { nr = 1./nr; }
    fAnalysisManager->GetH1(fHisto[0])->scale(nr);

    nr = sum*fStepR;
    if(nr > 0.0) { nr = 1./nr; }
    fAnalysisManager->GetH1(fHisto[1])->scale(nr);

    fAnalysisManager->GetH1(fHisto[3])
      ->scale(1000.0*cm3/(pi*fAbsorberR*fAbsorberR*fStepZ));
    fAnalysisManager->GetH1(fHisto[4])
      ->scale(1000.0*cm3*fVolumeR[0]/fStepZ);

    // Write histogram file
    if(!fAnalysisManager->Write()) {
      G4cout   << "Histo::Save: FATAL ERROR writing ROOT file" << G4endl;
      exit(1);
    }
    G4cout << "### Histo::Save: Histograms are saved" << G4endl;
    if(fAnalysisManager->CloseFile() && fVerbose > 0) {
      G4cout << "                 File is closed" << G4endl;
    }
    delete fAnalysisManager;
    fAnalysisManager = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::BookHisto()
{
  G4String nam = fHistName + ".root";

  if(!fAnalysisManager->OpenFile(nam)) {
    G4cout << "Histo::bookHisto: ERROR open file <" << nam << ">" << G4endl;
    exit(0); 
  }
  G4cout << "### Histo::Save: Opended file <" << nam << ">  for " 
         << fNHisto << " histograms " << G4endl;

  // Creating an 1-dimensional histograms in the root directory of the tree
  fHisto[0] = fAnalysisManager->CreateH1("10",
    "Energy deposit at radius (mm) normalised on 1st channel",fNBinsR,0.,fAbsorberR/mm);

  fHisto[1] = fAnalysisManager->CreateH1("11",
    "Energy deposit at radius (mm) normalised to integral",fNBinsR,0.,fAbsorberR/mm);

  fHisto[2] = fAnalysisManager->CreateH1("12",
    "Energy deposit (MeV/kg/electron) at radius (mm)",fNBinsR,0.,fAbsorberR/mm);

  fHisto[3] = fAnalysisManager->CreateH1("13",
    "Energy profile (MeV/kg/electron) over Z (mm)",fNBinsZ,0.,fAbsorberZ/mm);

  fHisto[4] = fAnalysisManager->CreateH1("14",
    "Energy profile (MeV/kg/electron) over Z (mm) at Central Voxel",fNBinsZ,0.,fAbsorberZ/mm);

  fHisto[5] = fAnalysisManager->CreateH1("15",
    "Energy (MeV) of fGamma produced in the target",fNBinsE,0.,fMaxEnergy/MeV);

  fHisto[6] = fAnalysisManager->CreateH1("16",
    "Energy (MeV) of fGamma before phantom",fNBinsE,0.,fMaxEnergy/MeV);

  fHisto[7] = fAnalysisManager->CreateH1("17",
    "Energy (MeV) of electrons produced in phantom",fNBinsE,0.,fMaxEnergy/MeV);

  fHisto[8] = fAnalysisManager->CreateH1("18",
    "Energy (MeV) of electrons produced in target",fNBinsE,0.,fMaxEnergy/MeV);

  fHisto[9] = fAnalysisManager->CreateH1("19",
    "Gamma Energy Fluence (MeV/cm2) at radius(mm) in front of phantom",fNBinsR,0.,fAbsorberR/mm);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddTargetPhoton(const G4DynamicParticle* p)
{
  ++fNgamTar; 
  if(fAnalysisManager) { 
    fAnalysisManager->FillH1(fHisto[5],p->GetKineticEnergy()/MeV,1.0); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddPhantomPhoton(const G4DynamicParticle*)
{
  ++fNgamPh; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddTargetElectron(const G4DynamicParticle* p)
{
  ++fNeTar; 
  if(fAnalysisManager) { 
    fAnalysisManager->FillH1(fHisto[8],p->GetKineticEnergy()/MeV,1.0); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddPhantomElectron(const G4DynamicParticle* p)
{
  ++fNePh;
  if(fAnalysisManager) { 
    fAnalysisManager->FillH1(fHisto[7],p->GetKineticEnergy()/MeV,1.0); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::ScoreNewTrack(const G4Track* aTrack)
{
  //Save primary parameters
  const G4ParticleDefinition* particle = aTrack->GetParticleDefinition();
  G4int pid = aTrack->GetParentID();
  G4VPhysicalVolume* pv = aTrack->GetVolume();
  const G4DynamicParticle* dp = aTrack->GetDynamicParticle();

  //primary particle
  if(0 == pid) {

    ++fNevt; 
    if(1 < fVerbose) {
      G4ThreeVector pos = aTrack->GetVertexPosition();
      G4ThreeVector dir = aTrack->GetMomentumDirection();
      G4cout << "TrackingAction: Primary "
             << particle->GetParticleName()
             << " Ekin(MeV)= " 
             << aTrack->GetKineticEnergy()/MeV
             << "; pos= " << pos << ";  dir= " << dir << G4endl;
    }

    // secondary electron
  } else if (0 < pid && particle == fElectron) {
    if(1 < fVerbose) {
      G4cout << "TrackingAction: Secondary electron " << G4endl;
    }
    AddElectron(dp);
    if(pv == fPhantom)                        { AddPhantomElectron(dp); }
    else if(pv == fTarget1 || pv == fTarget2) { AddTargetElectron(dp); }

    // secondary positron
  } else if (0 < pid && particle == fPositron) {
    if(1 < fVerbose) {
      G4cout << "TrackingAction: Secondary positron " << G4endl;
    }
    AddPositron(dp);

    // secondary gamma
  } else if (0 < pid && particle == fGamma) {
    if(1 < fVerbose) {
      G4cout << "TrackingAction: Secondary gamma; parentID= " << pid
             << " E= " << aTrack->GetKineticEnergy() << G4endl;
    }
    AddPhoton(dp);
    if(pv == fPhantom)                        { AddPhantomPhoton(dp); }
    else if(pv == fTarget1 || pv == fTarget2) { AddTargetPhoton(dp); }

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddPhantomGamma(G4double e, G4double r)
{
  e /= MeV;
  fSumR += e;
  G4int bin = (G4int)(e/fStepE);
  if(bin >= fNBinsE) { bin = fNBinsE-1; }
  fGammaE[bin] += e;
  G4int bin1 = (G4int)(r/fStepR);
  if(bin1 >= fNBinsR) { bin1 = fNBinsR-1; }
  if(fAnalysisManager) {
    fAnalysisManager->FillH1(fHisto[6],e,1.0);
    fAnalysisManager->FillH1(fHisto[9],r,e*fVolumeR[bin1]);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddPhantomStep(G4double edep, G4double r1, G4double z1, 
                           G4double r2, G4double z2,
                           G4double r0, G4double z0)
{
  ++fNstep;
  G4int nzbin = (G4int)(z0/fStepZ);
  if(fVerbose > 1) {
    G4cout << "Histo: edep(MeV)= " << edep/MeV << " at binz= " << nzbin
           << " r1= " << r1 << " z1= " << z1
           << " r2= " << r2 << " z2= " << z2
           << " r0= " << r0 << " z0= " << z0
           << G4endl;
  }
  if(nzbin == fScoreBin) {
    G4int bin  = (G4int)(r0/fStepR);
    if(bin >= fNBinsR) { bin = fNBinsR-1; }
    double w = edep*fVolumeR[bin];
    fEdep[bin] += w;
    if(fAnalysisManager) {
      fAnalysisManager->FillH1(fHisto[0],r0,w); 
      fAnalysisManager->FillH1(fHisto[1],r0,w); 
      fAnalysisManager->FillH1(fHisto[2],r0,w); 
    }
  }
  G4int bin1 = (G4int)(z1/fStepZ);
  if(bin1 >= fNBinsZ) { bin1 = fNBinsZ-1; }
  G4int bin2 = (G4int)(z2/fStepZ);
  if(bin2 >= fNBinsZ) { bin2 = fNBinsZ-1; }
  if(bin1 == bin2) {
    if(fAnalysisManager) {
      fAnalysisManager->FillH1(fHisto[3],z0,edep); 
      if(r1 < fStepR) {
        G4double w = edep;
        if(r2 > fStepR) { w *= (fStepR - r1)/(r2 - r1); }
        fAnalysisManager->FillH1(fHisto[4],z0,w);
      }
    }
  } else {
    G4int bin;

    if(bin2 < bin1) {
      bin = bin2;
      G4double z = z2;
      bin2 = bin1;
      z2 = z1;
      bin1 = bin;
      z1 = z;
    }
    G4double zz1 = z1;
    G4double zz2 = (bin1+1)*fStepZ;
    G4double rr1 = r1;
    G4double dz  = z2 - z1;
    G4double dr  = r2 - r1;
    G4double rr2 = r1 + dr*(zz2-zz1)/dz;
    for(bin=bin1; bin<=bin2; bin++) {
      if(fAnalysisManager) {
        G4double de = edep*(zz2 - zz1)/dz;
        G4double zf = (zz1+zz2)*0.5;
        { fAnalysisManager->FillH1(fHisto[3],zf,de); }
        if(rr1 < fStepR) {
          G4double w = de;
          if(rr2 > fStepR) w *= (fStepR - rr1)/(rr2 - rr1);
          { fAnalysisManager->FillH1(fHisto[4],zf,w); }
        }
      }
      zz1 = zz2;
      zz2 = std::min(z2, zz1+fStepZ);
      rr1 = rr2;
      rr2 = rr1 + dr*(zz2 - zz1)/dz;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

