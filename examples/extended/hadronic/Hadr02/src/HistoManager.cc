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
// $Id: HistoManager.cc,v 1.17 2010-01-13 15:53:44 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
//
// Author:      V.Ivanchenko 30/01/01
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
// 16.11.2006 Add beamFlag (V.Ivanchenko)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HistoManager.hh"
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
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "IonDPMJETPhysics.hh"
#include "G4IonFTFPBinaryCascadePhysics.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4BuilderType.hh"
 
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
{
  verbose=  0;
  nSlices   = 100;
  nBinsE    = 100;
  nHisto    = 20;
  length    = 300.*mm;
  edepMax   = 1.0*GeV;
  beamFlag  = true;
  material  = 0;
  elm       = 0;
  histo     = new Histo(verbose);
  neutron   = G4Neutron::Neutron();
  physList  = 0;
  ionPhysics= 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete histo;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::bookHisto()
{
  histo->add1D("0","Energy deposition (MeV/mm/event) in the target",
	       nSlices,0.0,length/mm,MeV/mm);
  histo->add1D("1","Log10 Energy (GeV) of gammas",nBinsE,-5.,5.,1.0);
  histo->add1D("2","Log10 Energy (GeV) of electrons",nBinsE,-5.,5.,1.0);
  histo->add1D("3","Log10 Energy (GeV) of positrons",nBinsE,-5.,5.,1.0);
  histo->add1D("4","Log10 Energy (GeV) of protons",nBinsE,-5.,5.,1.0);
  histo->add1D("5","Log10 Energy (GeV) of neutrons",nBinsE,-5.,5.,1.0);
  histo->add1D("6","Log10 Energy (GeV) of charged pions",nBinsE,-4.,6.,1.0);
  histo->add1D("7","Log10 Energy (GeV) of pi0",nBinsE,-4.,6.,1.0);
  histo->add1D("8","Log10 Energy (GeV) of charged kaons",nBinsE,-4.,6.,1.0);
  histo->add1D("9","Log10 Energy (GeV) of neutral kaons",nBinsE,-4.,6.,1.0);
  histo->add1D("10","Log10 Energy (GeV) of deuterons and tritons",nBinsE,-5.,5.,1.0);
  histo->add1D("11","Log10 Energy (GeV) of He3 and alpha",nBinsE,-5.,5.,1.0);
  histo->add1D("12","Log10 Energy (GeV) of Generic Ions",nBinsE,-5.,5.,1.0);
  histo->add1D("13","Log10 Energy (GeV) of muons",nBinsE,-4.,6.,1.0);
  histo->add1D("14","Log10 Energy (GeV) of pi+",nBinsE,-4.,6.,1.0);
  histo->add1D("15","Log10 Energy (GeV) of pi-",nBinsE,-4.,6.,1.0);
  histo->add1D("16","X Section (mb) of Secondary Fragments Z with E>1 GeV (mb)",25,0.5,25.5,1.0);
  histo->add1D("17","Secondary Fragment A E>1 GeV",50,0.5,50.5,1.0);
  histo->add1D("18","Secondary Fragment Z E<1 GeV",25,0.5,25.5,1.0);
  histo->add1D("19","Secondary Fragment A E<1 GeV",50,0.5,50.5,1.0);
  histo->add1D("20","X Section (mb) of Secondary Fragments Z (mb) ",25,0.5,25.5,1.0);
  histo->add1D("21","Secondary Fragment A ",50,0.5,50.5,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  absZ0       = -0.5*length;
  n_evt       = 0;
  n_elec      = 0;
  n_posit     = 0;
  n_gam       = 0;
  n_step      = 0;
  n_ions      = 0;
  n_deut      = 0;
  n_alpha     = 0;
  n_kaons     = 0;
  n_muons     = 0; 
  n_cpions    = 0;
  n_pi0       = 0;
  n_neutron   = 0;
  n_proton    = 0;
  n_aproton   = 0;

  edepSum     = 0.0;
  edepSum2    = 0.0;

  bookHisto();
  histo->book();

  if(verbose > 0) {
    G4cout << "HistoManager: Histograms are booked and run has been started"
           <<G4endl<<"  BeginOfRun (After histo->book)"<< G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun()
{

  G4cout << "HistoManager: End of run actions are started" << G4endl;

  // Average values
  G4cout<<"========================================================"<<G4endl;

  G4double x = (G4double)n_evt;
  if(n_evt > 0) { x = 1.0/x; }

  G4double xe = x*(G4double)n_elec;
  G4double xg = x*(G4double)n_gam;
  G4double xp = x*(G4double)n_posit;
  G4double xs = x*(G4double)n_step;
  G4double xn = x*(G4double)n_neutron;
  G4double xpn = x*(G4double)n_proton;
  G4double xap = x*(G4double)n_aproton;
  G4double xpc = x*(G4double)n_cpions;
  G4double xp0 = x*(G4double)n_pi0;
  G4double xpk = x*(G4double)n_kaons;
  G4double xpm = x*(G4double)n_muons;
  G4double xid = x*(G4double)n_deut;
  G4double xia = x*(G4double)n_alpha;
  G4double xio = x*(G4double)n_ions;

  edepSum  *= x;
  edepSum2 *= x;
  edepSum2 -= edepSum*edepSum;
  if(edepSum2 > 0.0) edepSum2 = std::sqrt(edepSum2);
  else               edepSum2 = 0.0;

  G4cout                         << "Beam particle                        "
				 << primaryDef->GetParticleName() <<G4endl;
  G4cout                         << "Beam Energy(GeV)                     " 
				 << primaryKineticEnergy/GeV <<G4endl;
  G4cout                         << "Number of events                     " << n_evt <<G4endl;
  G4cout << std::setprecision(4) << "Average energy deposit (GeV)         " << edepSum/GeV 
	 << "   RMS(GeV) " << edepSum2/GeV << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps              " << xs << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma              " << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e-                 " << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+                 " << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of neutrons           " << xn << G4endl;
  G4cout << std::setprecision(4) << "Average number of protons            " << xpn << G4endl;
  G4cout << std::setprecision(4) << "Average number of antiprotons        " << xap << G4endl;
  G4cout << std::setprecision(4) << "Average number of pi+ & pi-          " << xpc << G4endl;
  G4cout << std::setprecision(4) << "Average number of pi0                " << xp0 << G4endl;
  G4cout << std::setprecision(4) << "Average number of kaons              " << xpk << G4endl;
  G4cout << std::setprecision(4) << "Average number of muons              " << xpm << G4endl;
  G4cout << std::setprecision(4) << "Average number of deuterons+tritons  " << xid << G4endl;
  G4cout << std::setprecision(4) << "Average number of He3+alpha          " << xia << G4endl;
  G4cout << std::setprecision(4) << "Average number of ions               " << xio << G4endl;
  G4cout<<"========================================================"<<G4endl;
  G4cout<<G4endl;

  // normalise histograms
  for(G4int i=0; i<nHisto; i++) {
    histo->scale(i,x);
  }
  // will work only for pure material - 1 element
  //  G4cout << material << G4endl; 
  G4double F= 1000/(length*barn*material->GetTotNbOfAtomsPerVolume());
  if(F > 0.0) { 
    histo->scale(16, F); 
    histo->scale(20, F); 
  }

  if(verbose > 1) histo->print(0);
  histo->save();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfEvent()
{
  ++n_evt;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfEvent()
{
  edepSum  += edepEvt;
  edepSum2 += edepEvt*edepEvt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::ScoreNewTrack(const G4Track* track)
{
  const G4ParticleDefinition* pd = track->GetDefinition();
  currentDef = pd;
  G4String name = pd->GetParticleName();
  G4double e = track->GetKineticEnergy();
  currentKinEnergy = e;

  // Primary track
  if(0 == track->GetParentID()) {

    primaryKineticEnergy = e;
    primaryDef = pd;
    G4ThreeVector dir = track->GetMomentumDirection();
    if(1 < verbose) 
      G4cout << "### Primary " << name 
	     << " kinE(GeV)= " << e/GeV
	     << "; m(GeV)= " << pd->GetPDGMass()/GeV
	     << "; pos(mm)= " << track->GetPosition()/mm 
	     << ";  dir= " << track->GetMomentumDirection() 
	     << G4endl;
 
    // Secondary track
  } else {
    if(1 < verbose) 
      G4cout << "=== Secondary " << name 
	     << " kinE(GeV)= " << e/GeV
	     << "; m(GeV)= " << pd->GetPDGMass()/GeV
	     << "; pos(mm)= " << track->GetPosition()/mm 
	     << ";  dir= " << track->GetMomentumDirection() 
	     << G4endl;
    e = std::log10(e/GeV);
    if(pd == G4Gamma::Gamma()) {
      n_gam++;
      histo->fill(1,e,1.0);
    } else if ( pd == G4Electron::Electron()) {
      n_elec++;
      histo->fill(2,e,1.0);
    } else if ( pd == G4Positron::Positron()) {
      n_posit++;
      histo->fill(3,e,1.0);
    } else if ( pd == G4Proton::Proton()) {
      n_proton++;
      histo->fill(4,e,1.0);
    } else if ( pd == neutron) {
      n_neutron++;
      histo->fill(5,e,1.0);
    } else if ( pd == G4AntiProton::AntiProton()) {
      n_aproton++;
    } else if ( pd == G4PionPlus::PionPlus() ) {
      n_cpions++;
      histo->fill(6,e,1.0);
      histo->fill(14,e,1.0);

    } else if ( pd == G4PionMinus::PionMinus()) {
      n_cpions++;
      histo->fill(6,e,1.0);
      histo->fill(15,e,1.0);

    } else if ( pd == G4PionZero::PionZero()) {
      n_pi0++;
      histo->fill(7,e,1.0);
    } else if ( pd == G4KaonPlus::KaonPlus() || pd == G4KaonMinus::KaonMinus()) {
      n_kaons++;
      histo->fill(8,e,1.0);
    } else if ( pd == G4KaonZeroShort::KaonZeroShort() || pd == G4KaonZeroLong::KaonZeroLong()) {
      n_kaons++;
      histo->fill(9,e,1.0);
    } else if ( pd == G4Deuteron::Deuteron() || pd == G4Triton::Triton()) {
      n_deut++;
      histo->fill(10,e,1.0);
    } else if ( pd == G4He3::He3() || pd == G4Alpha::Alpha()) {
      n_alpha++;
      histo->fill(11,e,1.0);
    } else if ( pd->GetParticleType() == "nucleus") {
      n_ions++;
      histo->fill(12,e,1.0);
      G4double Z = pd->GetPDGCharge()/eplus;
      G4double A = (G4double)pd->GetBaryonNumber();
      if(e > 0.0) {
	histo->fill(16,Z,1.0);
	histo->fill(17,A,1.0);
      } else {
	histo->fill(18,Z,1.0);
	histo->fill(19,A,1.0);
      }
      histo->fill(20,Z,1.0);
      histo->fill(21,A,1.0);
    } else if ( pd == G4MuonPlus::MuonPlus() || pd == G4MuonMinus::MuonMinus()) {
      n_muons++;
      histo->fill(13,e,1.0);    
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddTargetStep(const G4Step* step)
{
  ++n_step;
  G4double edep = step->GetTotalEnergyDeposit();
  if(edep > 0.0) { 
    const G4Track* track = step->GetTrack();
    currentDef = track->GetDefinition(); 
    currentKinEnergy = track->GetKineticEnergy();

    G4ThreeVector pos = 
      (step->GetPreStepPoint()->GetPosition() +
       step->GetPostStepPoint()->GetPosition())*0.5;

    G4double z = pos.z() - absZ0;

    // scoring
    edepEvt += edep;
    histo->fill(0,z,edep);
 
    if(1 < verbose) {
      G4cout << "HistoManager::AddEnergy: e(keV)= " << edep/keV
	     << "; z(mm)= " << z/mm
	     << "; step(mm)= " << step->GetStepLength()/mm
	     << " by " << currentDef->GetParticleName()
	     << " E(GeV)= " << currentKinEnergy/GeV
	     << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetVerbose(G4int val)        
{
  verbose = val; 
  histo->setVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetTargetMaterial(const G4Material* mat)         
{
  if(mat) {
    material = mat;
    elm = (*(material->GetElementVector()))[0];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::Fill(G4int id, G4double x, G4double w)
{
  histo->fill(id, x, w);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetIonPhysics(const G4String& nam)
{
  if(ionPhysics) {
    G4cout << "### HistoManager WARNING: Ion Physics is already defined: <"
	   << nam << "> is ignored!" << G4endl;
  } else if(nam == "DPMJET") {
    ionPhysics = new IonDPMJETPhysics(false);
    physList->ReplacePhysics(ionPhysics);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "### SetIonPhysics: Ion Physics DPMJET/Binary is added"
	   << G4endl;
  } else if(nam == "FTF") {
    ionPhysics = new G4IonFTFPBinaryCascadePhysics();
    physList->ReplacePhysics(ionPhysics);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "### SetIonPhysics: Ion Physics FTFP/Binary is added"
	   << G4endl;
  } else {
    G4cout << "### HistoManager WARNING: Ion Physics <"
	   << nam << "> is unknown!" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

