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
/// \file hadronic/Hadr01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
//
// Author:      V.Ivanchenko 30/01/01
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
// 16.11.2006 Add beamFlag (V.Ivanchenko)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HistoManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZeroLong.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "Histo.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::GetPointer()
{
  if(!fManager) {
    static HistoManager manager;
    fManager = &manager;
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::HistoManager()
: fPrimaryDef(nullptr),
  fEdepMax(1.0*GeV),
  fLength (300.*mm),
  fPrimaryKineticEnergy(0.0),  
  fVerbose(0),  
  fNBinsE (100),
  fNSlices(300),
  fNHisto (28),
  fBeamFlag(true),
  fHistoBooked(false)
{
  fHisto     = new Histo();
  fHisto->SetVerbose(fVerbose);
  BookHisto();
  fNeutron   = G4Neutron::Neutron();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete fHisto;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BookHisto()
{
  fHistoBooked = true;
  fHisto->Add1D("1","Energy deposition (MeV/mm/event) in the target",
               fNSlices,0.0,fLength/mm,MeV/mm);
  fHisto->Add1D("2","Log10 Energy (MeV) of gammas",fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("3","Log10 Energy (MeV) of electrons",fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("4","Log10 Energy (MeV) of positrons",fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("5","Log10 Energy (MeV) of protons",fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("6","Log10 Energy (MeV) of neutrons",fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("7","Log10 Energy (MeV) of charged pions",fNBinsE,-4.,6.,1.0);
  fHisto->Add1D("8","Log10 Energy (MeV) of pi0",fNBinsE,-4.,6.,1.0);
  fHisto->Add1D("9","Log10 Energy (MeV) of charged kaons",fNBinsE,-4.,6.,1.0);
  fHisto->Add1D("10","Log10 Energy (MeV) of neutral kaons",fNBinsE,-4.,6.,1.0);
  fHisto->Add1D("11","Log10 Energy (MeV) of deuterons and tritons",
                fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("12","Log10 Energy (MeV) of He3 and alpha",fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("13","Log10 Energy (MeV) of Generic Ions",fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("14","Log10 Energy (MeV) of muons",fNBinsE,-4.,6.,1.0);
  fHisto->Add1D("15","log10 Energy (MeV) of side-leaked neutrons",
                fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("16","log10 Energy (MeV) of forward-leaked neutrons",
                fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("17","log10 Energy (MeV) of backward-leaked neutrons",
                fNBinsE,-5.,5.,1.0);
  fHisto->Add1D("18","log10 Energy (MeV) of leaking protons",
                fNBinsE,-4.,6.,1.0);
  fHisto->Add1D("19","log10 Energy (MeV) of leaking charged pions",
                fNBinsE,-4.,6.,1.0);
  fHisto->Add1D("20","Log10 Energy (MeV) of pi+",fNBinsE,-4.,6.,1.0);
  fHisto->Add1D("21","Log10 Energy (MeV) of pi-",fNBinsE,-4.,6.,1.0);
  fHisto->Add1D("22",
                "Energy deposition in the target normalized to beam energy",
                110,0.0,1.1,1.0);
  fHisto->Add1D("23",
                "EM energy deposition in the target normalized to beam energy",
                110,0.0,1.1,1.0);
  fHisto->Add1D("24",
               "Pion energy deposition in the target normalized to beam energy",
                110,0.0,1.1,1.0);
  fHisto->Add1D("25",
             "Proton energy deposition in the target normalized to beam energy",
                110,0.0,1.1,1.0);
  fHisto->Add1D("26","Energy (MeV) of pi+",fNBinsE,0.,500.,1.0);
  fHisto->Add1D("27","Energy (MeV) of pi-",fNBinsE,0.,500.,1.0);
  fHisto->Add1D("28","Energy (MeV) of pi0",fNBinsE,0.,500.,1.0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  fAbsZ0       = -0.5*fLength;
  fNevt       = 0;
  fNelec      = 0;
  fNposit     = 0;
  fNgam       = 0;
  fNstep      = 0;
  fNprot_leak = 0;
  fNpiofNleak = 0;
  fNions      = 0;
  fNdeut      = 0;
  fNalpha     = 0;
  fNkaons     = 0;
  fNmuons     = 0;
  fNcpions    = 0;
  fNpi0       = 0;
  fNneutron   = 0;
  fNproton    = 0;
  fNaproton   = 0;
  fNneu_forw  = 0;
  fNneu_leak  = 0;
  fNneu_back  = 0;

  fEdepSum     = 0.0;
  fEdepSum2    = 0.0;

  if(!fHistoBooked) { BookHisto(); }
  fHisto->Book();

  if(fVerbose > 0) {
    G4cout << "HistoManager: Histograms are booked and run has been started"
           <<G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun()
{

  G4cout << "HistoManager: End of run actions are started" << G4endl;

  // Average values
  G4cout<<"========================================================"<<G4endl;

  G4double x = (G4double)fNevt;
  if(fNevt > 0) { x = 1.0/x; }

  G4double xe = x*(G4double)fNelec;
  G4double xg = x*(G4double)fNgam;
  G4double xp = x*(G4double)fNposit;
  G4double xs = x*(G4double)fNstep;
  G4double xn = x*(G4double)fNneutron;
  G4double xpn = x*(G4double)fNproton;
  G4double xap = x*(G4double)fNaproton;
  G4double xnf = x*(G4double)fNneu_forw;
  G4double xnb = x*(G4double)fNneu_leak;
  G4double xnbw= x*(G4double)fNneu_back;
  G4double xpl = x*(G4double)fNprot_leak;
  G4double xal = x*(G4double)fNpiofNleak;
  G4double xpc = x*(G4double)fNcpions;
  G4double xp0 = x*(G4double)fNpi0;
  G4double xpk = x*(G4double)fNkaons;
  G4double xpm = x*(G4double)fNmuons;
  G4double xid = x*(G4double)fNdeut;
  G4double xia = x*(G4double)fNalpha;
  G4double xio = x*(G4double)fNions;

  fEdepSum  *= x;
  fEdepSum2 *= x;
  fEdepSum2 -= fEdepSum*fEdepSum;
  if(fEdepSum2 > 0.0) { fEdepSum2 = std::sqrt(fEdepSum2); }
  else                { fEdepSum2 = 0.0; }

  G4cout                         << "Beam particle                        "
                                 << fPrimaryDef->GetParticleName() <<G4endl;
  G4cout                         << "Beam Energy(MeV)                     " 
                                 << fPrimaryKineticEnergy/MeV <<G4endl;
  G4cout                         << "Number of events                     " 
                                 << fNevt <<G4endl;
  G4cout << std::setprecision(4) << "Average energy deposit (MeV)         " 
         << fEdepSum/MeV 
         << "   RMS(MeV) " << fEdepSum2/MeV << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps              " 
         << xs << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma              " 
         << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e-                 " 
         << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+                 " 
         << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of neutrons           " 
         << xn << G4endl;
  G4cout << std::setprecision(4) << "Average number of protons            " 
         << xpn << G4endl;
  G4cout << std::setprecision(4) << "Average number of antiprotons        " 
         << xap << G4endl;
  G4cout << std::setprecision(4) << "Average number of pi+ & pi-          " 
         << xpc << G4endl;
  G4cout << std::setprecision(4) << "Average number of pi0                " 
         << xp0 << G4endl;
  G4cout << std::setprecision(4) << "Average number of kaons              " 
         << xpk << G4endl;
  G4cout << std::setprecision(4) << "Average number of muons              " 
         << xpm << G4endl;
  G4cout << std::setprecision(4) << "Average number of deuterons+tritons  " 
         << xid << G4endl;
  G4cout << std::setprecision(4) << "Average number of He3+alpha          " 
         << xia << G4endl;
  G4cout << std::setprecision(4) << "Average number of ions               " 
         << xio << G4endl;
  G4cout << std::setprecision(4) << "Average number of forward neutrons   " 
         << xnf << G4endl;
  G4cout << std::setprecision(4) << "Average number of reflected neutrons " 
         << xnb << G4endl;
  G4cout << std::setprecision(4) << "Average number of leaked neutrons    " 
         << xnbw << G4endl;
  G4cout << std::setprecision(4) << "Average number of proton leak        " 
         << xpl << G4endl;
  G4cout << std::setprecision(4) << "Average number of pion leak          " 
         << xal << G4endl;
  G4cout<<"========================================================"<<G4endl;
  G4cout<<G4endl;

  // normalise histograms
  for(G4int i=0; i<fNHisto; i++) { 
    fHisto->ScaleH1(i,x);
  }

  fHisto->Save();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfEvent()
{
  fEdepEvt = 0.0;
  fEdepEM  = 0.0;
  fEdepPI  = 0.0;
  fEdepP   = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfEvent()
{
  fEdepSum  += fEdepEvt;
  fEdepSum2 += fEdepEvt*fEdepEvt;
  fHisto->Fill(21,fEdepEvt/fPrimaryKineticEnergy,1.0);
  fHisto->Fill(22,fEdepEM/fPrimaryKineticEnergy,1.0);
  fHisto->Fill(23,fEdepPI/fPrimaryKineticEnergy,1.0);
  fHisto->Fill(24,fEdepP/fPrimaryKineticEnergy,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::ScoreNewTrack(const G4Track* track)
{
  const G4ParticleDefinition* pd = track->GetDefinition();
  G4String name = pd->GetParticleName();
  G4double ee = track->GetKineticEnergy();
  G4double e = ee;

  // Primary track
  if(0 == track->GetParentID()) {

    fNevt++;
    fPrimaryKineticEnergy = e;
    fPrimaryDef = pd;
    G4ThreeVector dir = track->GetMomentumDirection();
    if(1 < fVerbose) 
      G4cout << "### Primary " << name 
             << " kinE(MeV)= " << e/MeV
             << "; m(MeV)= " << pd->GetPDGMass()/MeV
             << "; pos(mm)= " << track->GetPosition()/mm 
             << ";  dir= " << track->GetMomentumDirection() 
             << G4endl;

    // Secondary track
  } else {
    if(1 < fVerbose) 
      G4cout << "=== Secondary " << name 
             << " kinE(MeV)= " << e/MeV
             << "; m(MeV)= " << pd->GetPDGMass()/MeV
             << "; pos(mm)= " << track->GetPosition()/mm 
             << ";  dir= " << track->GetMomentumDirection() 
             << G4endl;
    e = std::log10(e/MeV);
    if(pd == G4Gamma::Gamma()) {
      fNgam++;
      fHisto->Fill(1,e,1.0);
    } else if ( pd == G4Electron::Electron()) {
      fNelec++;
      fHisto->Fill(2,e,1.0);
    } else if ( pd == G4Positron::Positron()) {
      fNposit++;
      fHisto->Fill(3,e,1.0);
    } else if ( pd == G4Proton::Proton()) {
      fNproton++;
      fHisto->Fill(4,e,1.0);
    } else if ( pd == fNeutron) {
      fNneutron++;
      fHisto->Fill(5,e,1.0);
    } else if ( pd == G4AntiProton::AntiProton()) {
      fNaproton++;
    } else if ( pd == G4PionPlus::PionPlus() ) {
      fNcpions++;
      fHisto->Fill(6,e,1.0);
      fHisto->Fill(19,e,1.0);
      fHisto->Fill(25,ee,1.0);

    } else if ( pd == G4PionMinus::PionMinus()) {
      fNcpions++;
      fHisto->Fill(6,e,1.0);
      fHisto->Fill(20,e,1.0);
      fHisto->Fill(26,ee,1.0);

    } else if ( pd == G4PionZero::PionZero()) {
      fNpi0++;
      fHisto->Fill(7,e,1.0);
      fHisto->Fill(27,ee,1.0);

    } else if ( pd == G4KaonPlus::KaonPlus() || 
                pd == G4KaonMinus::KaonMinus()) {
      fNkaons++;
      fHisto->Fill(8,e,1.0);
    } else if ( pd == G4KaonZeroShort::KaonZeroShort() || 
                pd == G4KaonZeroLong::KaonZeroLong()) {
      fNkaons++;
      fHisto->Fill(9,e,1.0);
    } else if ( pd == G4Deuteron::Deuteron() || pd == G4Triton::Triton()) {
      fNdeut++;
      fHisto->Fill(10,e,1.0);
    } else if ( pd == G4He3::He3() || pd == G4Alpha::Alpha()) {
      fNalpha++;
      fHisto->Fill(11,e,1.0);
    } else if ( pd->GetParticleType() == "nucleus") {
      fNions++;
      fHisto->Fill(12,e,1.0);
    } else if ( pd == G4MuonPlus::MuonPlus() || 
                pd == G4MuonMinus::MuonMinus()) {
      fNmuons++;
      fHisto->Fill(13,e,1.0);    
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddTargetStep(const G4Step* step)
{
  fNstep++;
  G4double fEdep = step->GetTotalEnergyDeposit();
  if(1 < fVerbose) {
    G4cout << "TargetSD::ProcessHits: beta1= " 
           << step->GetPreStepPoint()->GetVelocity()/c_light
           << "  beta2= " << step->GetPostStepPoint()->GetVelocity()/c_light
           << " weight= " << step->GetTrack()->GetWeight() 
           << G4endl;
  }
  if(fEdep >= DBL_MIN) { 
    const G4Track* track = step->GetTrack();

    G4ThreeVector pos = 
      (step->GetPreStepPoint()->GetPosition() +
       step->GetPostStepPoint()->GetPosition())*0.5;

    G4double z = pos.z() - fAbsZ0;

    // scoring
    fEdepEvt += fEdep;
    fHisto->Fill(0,z,fEdep);
    const G4ParticleDefinition* pd = track->GetDefinition();

    if(pd == G4Gamma::Gamma() || pd == G4Electron::Electron() 
       || pd == G4Positron::Positron()) {
      fEdepEM += fEdep;
    } else if ( pd == G4PionPlus::PionPlus() || 
                pd == G4PionMinus::PionMinus()) {
      fEdepPI += fEdep;
    } else if ( pd == G4Proton::Proton() || 
                pd == G4AntiProton::AntiProton()) {
      fEdepP  += fEdep;
    }

    if(1 < fVerbose) {
      G4cout << "HistoManager::AddEnergy: e(keV)= " << fEdep/keV
             << "; z(mm)= " << z/mm
             << "; step(mm)= " << step->GetStepLength()/mm
             << " by " << pd->GetParticleName()
             << " E(MeV)= " << track->GetKineticEnergy()/MeV
             << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddLeakingParticle(const G4Track* track)
{
  const G4ParticleDefinition* pd = track->GetDefinition();
  const G4StepPoint* sp = track->GetStep()->GetPreStepPoint(); 
  G4double e = std::log10(sp->GetKineticEnergy()/MeV);

  G4ThreeVector pos = sp->GetPosition();
  G4ThreeVector dir = sp->GetMomentumDirection();
  G4double x = pos.x();
  G4double y = pos.y();
  G4double z = pos.z();
 
  G4bool isLeaking = false;

  // Forward 
  if(z > -fAbsZ0 && dir.z() > 0.0) {
    isLeaking = true;
    if(pd == fNeutron) {
      ++fNneu_forw;
      fHisto->Fill(15,e,1.0);
    } else isLeaking = true;

    // Backward
  } else if (z < fAbsZ0 && dir.z() < 0.0) {
    isLeaking = true;
    if(pd == fNeutron) {
      ++fNneu_back;
      fHisto->Fill(16,e,1.0);
    } else isLeaking = true;

    // Side
  } else if (std::abs(z) <= -fAbsZ0 && x*dir.x() + y*dir.y() > 0.0) {
    isLeaking = true;
    if(pd == fNeutron) {
      ++fNneu_leak;
      fHisto->Fill(14,e,1.0);
    } else isLeaking = true;
  }

  // protons and pions
  if(isLeaking) {
    if(pd == G4Proton::Proton()) {
      fHisto->Fill(17,e,1.0);
      ++fNprot_leak;
    } else if (pd == G4PionPlus::PionPlus() || 
               pd == G4PionMinus::PionMinus()) {
      fHisto->Fill(18,e,1.0);
      ++fNpiofNleak;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetVerbose(G4int val)        
{
  fVerbose = val; 
  fHisto->SetVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::Fill(G4int id, G4double x, G4double w)
{
  fHisto->Fill(id, x, w);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

