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
/// \file electromagnetic/TestEm9/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
//
// Author:      V.Ivanchenko 30/01/01
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HistoManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4EmProcessSubType.hh"
#include "G4VProcess.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4UnitsTable.hh"
#include "Histo.hh"
#include "EmAcceptance.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4SystemOfUnits.hh"
#include "G4GammaGeneralProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::fManager = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::GetPointer()
{
  if(nullptr == fManager) {
    fManager = new HistoManager();
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::HistoManager()
 : fGamma(G4Gamma::Gamma()),
   fElectron(G4Electron::Electron()),
   fPositron(G4Positron::Positron()),
   fHisto(new Histo())
{
  fVerbose = 1;
  fEvt1    = -1;
  fEvt2    = -1;
  fNmax    = 3;
  fMaxEnergy = 50.0*MeV;
  fBeamEnergy= 1.*GeV;
  fMaxEnergyAbs = 10.*MeV;
  fBinsE = 100;
  fBinsEA= 40;
  fBinsED= 100;
  fNHisto = 13;

  BookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete fHisto;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BookHisto()
{
  fHisto->Add1D("10","Evis/E0 in central crystal",fBinsED,0.0,1,1.0);
  fHisto->Add1D("11","Evis/E0 in 3x3",fBinsED,0.0,1.0,1.0);
  fHisto->Add1D("12","Evis/E0 in 5x5",fBinsED,0.0,1.0,1.0);
  fHisto->Add1D("13","Energy (MeV) of delta-electrons",
                fBinsE,0.0,fMaxEnergy,MeV);
  fHisto->Add1D("14","Energy (MeV) of gammas",fBinsE,0.0,fMaxEnergy,MeV);
  fHisto->Add1D("15","Energy (MeV) in abs1",fBinsEA,0.0,fMaxEnergyAbs,MeV);
  fHisto->Add1D("16","Energy (MeV) in abs2",fBinsEA,0.0,fMaxEnergyAbs,MeV);
  fHisto->Add1D("17","Energy (MeV) in abs3",fBinsEA,0.0,fMaxEnergyAbs,MeV);
  fHisto->Add1D("18","Energy (MeV) in abs4",fBinsEA,0.0,fMaxEnergyAbs,MeV);
  fHisto->Add1D("19","Number of vertex hits",20,-0.5,19.5,1.0);
  fHisto->Add1D("20","E1/E9 Ratio",fBinsED,0.0,1,1.0);
  fHisto->Add1D("21","E1/E25 Ratio",fBinsED,0.0,1.0,1.0);
  fHisto->Add1D("22","E9/E25 Ratio",fBinsED,0.0,1.0,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  // initilise scoring
  fEvt  = 0;
  fElec = 0;
  fPosit= 0;
  fGam  = 0;
  fStep = 0;
  fLowe = 0;

  for(G4int i=0; i<6; i++) {
    fStat[i] = 0;
    fEdep[i] = 0.0;
    fErms[i] = 0.0;
    if(i < 3) {
      fEdeptr[i] = 0.0;
      fErmstr[i] = 0.0;
    }
  }

  // initialise counters
  fBrem.resize(93,0.0);
  fPhot.resize(93,0.0);
  fComp.resize(93,0.0);
  fConv.resize(93,0.0);

  // initialise acceptance - by default is not applied
  for(G4int i=0; i<fNmax; i++) {
    fEdeptrue[i] = 1.0;
    fRmstrue[i]  = 1.0;
    fLimittrue[i]= 10.;
  }

  if(fHisto->IsActive()) { 
    for(G4int i=0; i<fNHisto; ++i) {fHisto->Activate(i, true); }
    fHisto->Book();

    if(fVerbose > 0) {
      G4cout << "HistoManager: Histograms are booked and run has been started"
             << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun(G4int runID)
{

  G4cout << "HistoManager: End of run actions are started   RunID# " 
         << runID << G4endl;
  G4String nam[6] = {"1x1", "3x3", "5x5", "E1/E9 ", "E1/E25", "E9/E25"};

  // average

  G4cout<<"================================================================="
        <<G4endl;
  G4double x = (G4double)fEvt;
  if(fEvt > 0) x = 1.0/x;
  G4int j;
  for(j=0; j<fNmax; j++) {

    // total mean
    fEdep[j] *= x;
    G4double y = fErms[j]*x - fEdep[j]*fEdep[j];
    if(y < 0.0) y = 0.0;
    fErms[j] = std::sqrt(y);

    // trancated mean
    G4double xx = G4double(fStat[j]);
    if(xx > 0.0) xx = 1.0/xx;
    fEdeptr[j] *= xx;
    y = fErmstr[j]*xx - fEdeptr[j]*fEdeptr[j];
    if(y < 0.0) y = 0.0;
    fErmstr[j] = std::sqrt(y);
  }
  G4double xe = x*(G4double)fElec;
  G4double xg = x*(G4double)fGam;
  G4double xp = x*(G4double)fPosit;
  G4double xs = x*fStep;

  G4double f = 100.*std::sqrt(fBeamEnergy/GeV);

  G4cout                         << "Number of events             "
                                 << fEvt <<G4endl;
  G4cout << std::setprecision(4) << "Average number of e-         "
         << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma      "
         << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+         "
         << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps      "
         << xs << G4endl;
  
  for(j=0; j<3; ++j) {
    G4double ex = fEdeptr[j];
    G4double sx = fErmstr[j];
    G4double xx= G4double(fStat[j]);
    if(xx > 0.0) xx = 1.0/xx;
    G4double r = sx*std::sqrt(xx);
    G4cout << std::setprecision(4) << "Edep " << nam[j]
           << " =                   " << ex
           << " +- " << r;
    if(ex > 0.1) G4cout << "  res=  " << f*sx/ex << " %   " << fStat[j];
    G4cout << G4endl;
  }
  if(fLimittrue[0] < 10. || fLimittrue[1] < 10. || fLimittrue[2] < 10.) {
    G4cout
      <<"===========  Mean values without trancating ====================="
      <<G4endl;
    for(j=0; j<fNmax; j++) {
      G4double ex = fEdep[j];
      G4double sx = fErms[j];
      G4double rx = sx*std::sqrt(x);
      G4cout << std::setprecision(4) << "Edep " << nam[j]
             << " =                   " << ex
             << " +- " << rx;
      if(ex > 0.0) G4cout << "  res=  " << f*sx/ex << " %";
      G4cout << G4endl;
    }
  }
  G4cout
    <<"===========  Ratios without trancating ==========================="<<G4endl;
  for(j=3; j<6; ++j) {
    G4double e = fEdep[j];
    G4double xx= G4double(fStat[j]);
    if(xx > 0.0) xx = 1.0/xx;
    e *= xx;
    G4double y = fErms[j]*xx - e*e;
    G4double r = 0.0;
    if(y > 0.0) r = std::sqrt(y*xx);
    G4cout << "  " << nam[j] << " =                   " << e
           << " +- " << r;
    G4cout << G4endl;
  }
  G4cout << std::setprecision(4) << "Beam Energy                  "
         << fBeamEnergy/GeV
         << " GeV" << G4endl;
  if(fLowe > 0)          G4cout << "Number of events E/E0<0.8    "
                                << fLowe << G4endl; 
  G4cout
    <<"=================================================================="
    <<G4endl;
  G4cout<<G4endl;

  // normalise histograms
  if(fHisto->IsActive()) { 
    for(G4int i=0; i<fNHisto; ++i) {
      fHisto->ScaleH1(i,x);
    }
    fHisto->Save();
  }
  if(0 < runID) { return; }

  // Acceptance only for the first run
  EmAcceptance acc;
  G4bool isStarted = false;
  for (j=0; j<fNmax; j++) {

    G4double ltrue = fLimittrue[j];
    if (ltrue < DBL_MAX) {
      if (!isStarted) {
        acc.BeginOfAcceptance("Crystal Calorimeter",fEvt);
        isStarted = true;
      }
      G4double etrue = fEdeptrue[j];
      G4double rtrue = fRmstrue[j];
      acc.EmAcceptanceGauss("Edep"+nam[j],fEvt,fEdeptr[j],etrue,rtrue,ltrue);
      acc.EmAcceptanceGauss("Erms"+nam[j],fEvt,fErmstr[j],rtrue,rtrue,
                            2.0*ltrue);
    }
  }
  if(isStarted) acc.EndOfAcceptance();

  // atom frequency
  G4cout << "   Z  bremsstrahlung photoeffect  compton    conversion" << G4endl;
  for(j=1; j<93; ++j) {
    G4int n1 = G4int(fBrem[j]*x);
    G4int n2 = G4int(fPhot[j]*x);
    G4int n3 = G4int(fComp[j]*x);
    G4int n4 = G4int(fConv[j]*x);
    if(n1 + n2 + n3 + n4 > 0) {
      G4cout << std::setw(4) << j << std::setw(12) << n1 
             << std::setw(12) << n2
             << std::setw(12) << n3 << std::setw(12) << n4 << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfEvent()
{
  ++fEvt;

  fEabs1  = 0.0;
  fEabs2  = 0.0;
  fEabs3  = 0.0;
  fEabs4  = 0.0;
  fEvertex.clear();
  fNvertex.clear();
  for (G4int i=0; i<25; i++) {
    fE[i] = 0.0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfEvent()
{
  G4double e9 = 0.0;
  G4double e25= 0.0;
  for (G4int i=0; i<25; i++) {
    fE[i] /= fBeamEnergy;
    e25 += fE[i];
    if( ( 6<=i &&  8>=i) || (11<=i && 13>=i) || (16<=i && 18>=i)) e9 += fE[i];
  }

  if(1 < fVerbose && e25 < 0.8) {
    ++fLowe;
    G4cout << "### in the event# " << fEvt << "  E25= " << e25 << G4endl;
  }

  // compute ratios
  G4double e0 = fE[12];
  G4double e19  = 0.0;
  G4double e125 = 0.0;
  G4double e925 = 0.0;
  if(e9 > 0.0) {
    e19 = e0/e9;
    e125 = e0/e25;
    e925 = e9/e25;
    fEdep[3] += e19;
    fErms[3] += e19*e19;
    fEdep[4] += e125;
    fErms[4] += e125*e125;
    fEdep[5] += e925;
    fErms[5] += e925*e925;
    fStat[3] += 1;
    fStat[4] += 1;
    fStat[5] += 1;
  }

  // Fill histo
  fHisto->Fill(0,e0,1.0);
  fHisto->Fill(1,e9,1.0);
  fHisto->Fill(2,e25,1.0);
  fHisto->Fill(5,fEabs1,1.0);
  fHisto->Fill(6,fEabs2,1.0);
  fHisto->Fill(7,fEabs3,1.0);
  fHisto->Fill(8,fEabs4,1.0);
  fHisto->Fill(9,G4double(fNvertex.size()),1.0);
  fHisto->Fill(10,e19,1.0);
  fHisto->Fill(11,e125,1.0);
  fHisto->Fill(12,e925,1.0);

  // compute sums
  fEdep[0] += e0;
  fErms[0] += e0*e0;
  fEdep[1] += e9;
  fErms[1] += e9*e9;
  fEdep[2] += e25;
  fErms[2] += e25*e25;

  // trancated mean
  if(std::abs(e0-fEdeptrue[0])<fRmstrue[0]*fLimittrue[0]) {
    fStat[0] += 1;
    fEdeptr[0] += e0;
    fErmstr[0] += e0*e0;
  }
  if(std::abs(e9-fEdeptrue[1])<fRmstrue[1]*fLimittrue[1]) {
    fStat[1] += 1;
    fEdeptr[1] += e9;
    fErmstr[1] += e9*e9;
  }
  if(std::abs(e25-fEdeptrue[2])<fRmstrue[2]*fLimittrue[2]) {
    fStat[2] += 1;
    fEdeptr[2] += e25;
    fErmstr[2] += e25*e25;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::ScoreNewTrack(const G4Track* aTrack)
{
  //Save primary parameters
  ResetTrackLength();
  const G4ParticleDefinition* particle = aTrack->GetDefinition();
  const G4DynamicParticle* dynParticle = aTrack->GetDynamicParticle();

  G4int pid = aTrack->GetParentID();
  G4double kinE = dynParticle->GetKineticEnergy();
  G4ThreeVector pos = aTrack->GetVertexPosition();

  // primary
  if(0 == pid) {

    G4double mass = 0.0;
    if(particle) {
      mass = particle->GetPDGMass();
    }

    G4ThreeVector dir = dynParticle->GetMomentumDirection();
    if(1 < fVerbose) {
      G4cout << "TrackingAction: Primary kinE(MeV)= " << kinE/MeV
           << "; m(MeV)= " << mass/MeV
           << "; pos= " << pos << ";  dir= " << dir << G4endl;
    }

    // secondary
  } else {
    const G4VProcess* proc = aTrack->GetCreatorProcess();
    G4int type = proc->GetProcessSubType();
   
    if(type == fBremsstrahlung) {
      auto elm = 
        static_cast<const G4VEnergyLossProcess*>(proc)->GetCurrentElement();
      if(nullptr != elm) {
        G4int Z = elm->GetZasInt();
        if(Z > 0 && Z < 93) { fBrem[Z] += 1.0; }
      }
    } else if(type == fPhotoElectricEffect) {
      auto elm = static_cast<const G4VEmProcess*>(proc)->GetCurrentElement();
      if(nullptr != elm) {
        G4int Z = elm->GetZasInt();
        if(Z > 0 && Z < 93) { fPhot[Z] += 1.0; }
      }
    } else if(type == fGammaConversion) {
      auto elm = static_cast<const G4VEmProcess*>(proc)->GetCurrentElement();
      if(nullptr != elm) {
        G4int Z = elm->GetZasInt();
        if(Z > 0 && Z < 93) { fConv[Z] += 1.0; }
      }
    } else if(type == fComptonScattering) {
      auto elm = static_cast<const G4VEmProcess*>(proc)->GetCurrentElement();
      if(nullptr != elm) {
        G4int Z = elm->GetZasInt();
        if(Z > 0 && Z < 93) { fComp[Z] += 1.0; }
      }
    }

    // delta-electron
    if (particle == fElectron) {
      if(1 < fVerbose) {
        G4cout << "TrackingAction: Secondary electron " << G4endl;
      }
      AddDeltaElectron(dynParticle);

    } else if (particle == fPositron) {
      if(1 < fVerbose) {
        G4cout << "TrackingAction: Secondary positron " << G4endl;
      }
      AddPositron(dynParticle);

    } else if (particle == fGamma) {
      if(1 < fVerbose) {
        G4cout << "TrackingAction: Secondary gamma; parentID= " << pid
               << " E= " << kinE << G4endl;
      }
      AddPhoton(dynParticle);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddEnergy(G4double edep, G4int volIndex, G4int copyNo)
{
  if(1 < fVerbose) {
    G4cout << "HistoManager::AddEnergy: e(keV)= " << edep/keV
           << "; volIdx= " << volIndex
           << "; copyNo= " << copyNo
           << G4endl;
  }
  if(0 == volIndex) {
    fE[copyNo] += edep;
  } else if (1 == volIndex) {
    fEabs1 += edep;
  } else if (2 == volIndex) {
    fEabs2 += edep;
  } else if (3 == volIndex) {
    fEabs3 += edep;
  } else if (4 == volIndex) {
    fEabs4 += edep;
  } else if (5 == volIndex) {
    G4int n = fNvertex.size();
    G4bool newPad = true;
    if (n > 0) {
      for(G4int i=0; i<n; i++) {
        if (copyNo == fNvertex[i]) {
          newPad = false;
          fEvertex[i] += edep;
          break;
        }
      }
    }
    if(newPad) {
      fNvertex.push_back(copyNo);
      fEvertex.push_back(edep);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddDeltaElectron(const G4DynamicParticle* elec)
{
  G4double e = elec->GetKineticEnergy()/MeV;
  if(e > 0.0) {
    ++fElec;
    fHisto->Fill(3,e,1.0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddPhoton(const G4DynamicParticle* ph)
{
  G4double e = ph->GetKineticEnergy()/MeV;
  if(e > 0.0) { 
    ++fGam;
    fHisto->Fill(4,e,1.0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetEdepAndRMS(G4int i, const G4ThreeVector& val)
{
  if(i<fNmax && i>=0) {
    if(val[0] > 0.0) fEdeptrue[i] = val[0];
    if(val[1] > 0.0) fRmstrue[i] = val[1];
    if(val[2] > 0.0) fLimittrue[i] = val[2];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
