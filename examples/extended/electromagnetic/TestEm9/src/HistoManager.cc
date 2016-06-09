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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::GetPointer()
{
  if(!fManager) {
    fManager = new HistoManager();
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::HistoManager()
{
  verbose = 1;
  nEvt1   = -1;
  nEvt2   = -1;
  nmax    = 3;
  histo   = new Histo();
  bookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete histo;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::bookHisto()
{
  maxEnergy = 50.0*MeV;
  beamEnergy= 1.*GeV;
  maxEnergyAbs = 10.*MeV;
  thKinE    = 1.*MeV;
  nBinsE = 100;
  nBinsEA= 40;
  nBinsED= 100;
  nTuple = false;
  nHisto = 13;

  // initialise acceptance
  for(G4int i=0; i<nmax; i++) {
    edeptrue[i] = 1.0;
    rmstrue[i]  = 1.0;
    limittrue[i]= DBL_MAX;
  }

  histo->add1D("10",
    "Evis/E0 in central crystal",nBinsED,0.0,1,1.0);

  histo->add1D("11",
    "Evis/E0 in 3x3",nBinsED,0.0,1.0,1.0);

  histo->add1D("12",
    "Evis/E0 in 5x5",nBinsED,0.0,1.0,1.0);

  histo->add1D("13",
    "Energy (MeV) of delta-electrons",nBinsE,0.0,maxEnergy,MeV);

  histo->add1D("14",
    "Energy (MeV) of gammas",nBinsE,0.0,maxEnergy,MeV);

  histo->add1D("15",
    "Energy (MeV) in abs1",nBinsEA,0.0,maxEnergyAbs,MeV);

  histo->add1D("16",
    "Energy (MeV) in abs2",nBinsEA,0.0,maxEnergyAbs,MeV);

  histo->add1D("17",
    "Energy (MeV) in abs3",nBinsEA,0.0,maxEnergyAbs,MeV);

  histo->add1D("18",
    "Energy (MeV) in abs4",nBinsEA,0.0,maxEnergyAbs,MeV);

  histo->add1D("19",
    "Number of vertex hits",20,-0.5,19.5,1.0);

  histo->add1D("20",
    "E1/E9 Ratio",nBinsED,0.0,1,1.0);

  histo->add1D("21",
    "E1/25 Ratio",nBinsED,0.0,1.0,1.0);

  histo->add1D("22",
    "E9/E25 Ratio",nBinsED,0.0,1.0,1.0);

  if(nTuple) {
    histo->addTuple( "100", "Dose deposite","float r, z, e" );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  // initilise scoring
  n_evt  = 0;
  n_elec = 0;
  n_posit= 0;
  n_gam  = 0;
  n_step = 0;
  n_lowe = 0;

  for(G4int i=0; i<6; i++) {
    stat[i] = 0;
    edep[i] = 0.0;
    erms[i] = 0.0;
    if(i < 3) {
      edeptr[i] = 0.0;
      ermstr[i] = 0.0;
    }
  }

  histo->book();

  brem.resize(93,0.0);
  phot.resize(93,0.0);
  comp.resize(93,0.0);
  conv.resize(93,0.0);

  if(verbose > 0) {
    G4cout << "HistoManager: Histograms are booked and run has been started"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun(G4int runID)
{

  G4cout << "HistoManager: End of run actions are started   RunID# " 
	 << runID << G4endl;
  G4String nam[6] = {"1x1", "3x3", "5x5", "E1/E9 ", "E1/E25", "E9/E25"};

  // average

  G4cout<<"================================================================="<<G4endl;
  G4double x = (G4double)n_evt;
  if(n_evt > 0) x = 1.0/x;
  G4int j;
  for(j=0; j<nmax; j++) {

    // total mean
    edep[j] *= x;
    G4double y = erms[j]*x - edep[j]*edep[j];
    if(y < 0.0) y = 0.0;
    erms[j] = std::sqrt(y);

    // trancated mean
    G4double xx = G4double(stat[j]);
    if(xx > 0.0) xx = 1.0/xx;
    edeptr[j] *= xx;
    y = ermstr[j]*xx - edeptr[j]*edeptr[j];
    if(y < 0.0) y = 0.0;
    ermstr[j] = std::sqrt(y);
  }
  G4double xe = x*(G4double)n_elec;
  G4double xg = x*(G4double)n_gam;
  G4double xp = x*(G4double)n_posit;
  G4double xs = x*n_step;

  G4double f = 100.*std::sqrt(beamEnergy/GeV);

  G4cout                         << "Number of events             " << n_evt <<G4endl;
  G4cout << std::setprecision(4) << "Average number of e-         " << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma      " << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+         " << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps      " << xs << G4endl;
  
  for(j=0; j<3; j++) {
    G4double e = edeptr[j];
    G4double s = ermstr[j];
    G4double xx= G4double(stat[j]);
    if(xx > 0.0) xx = 1.0/xx;
    G4double r = s*std::sqrt(xx);
    G4cout << std::setprecision(4) << "Edep " << nam[j] << " =                   " << e
           << " +- " << r;
    if(e > 0.0) G4cout << "  res=  " << f*s/e << " %";
    G4cout << G4endl;
  }
  if(limittrue[0] != DBL_MAX || limittrue[1] != DBL_MAX || limittrue[2] != DBL_MAX) {
    G4cout<<"===========  Mean values without trancating ====================="<<G4endl;
    for(j=0; j<nmax; j++) {
      G4double e = edep[j];
      G4double s = erms[j];
      G4double r = s*std::sqrt(x);
      G4cout << std::setprecision(4) << "Edep " << nam[j] << " =                   " << e
	     << " +- " << r;
      if(e > 0.0) G4cout << "  res=  " << f*s/e << " %";
      G4cout << G4endl;
    }
  }
  G4cout<<"===========  Ratios without trancating ==========================="<<G4endl;
  for(j=3; j<6; j++) {
    G4double e = edep[j];
    G4double xx= G4double(stat[j]);
    if(xx > 0.0) xx = 1.0/xx;
    e *= xx;
    G4double y = erms[j]*xx - e*e;
    G4double r = 0.0;
    if(y > 0.0) r = std::sqrt(y*xx);
    G4cout << "  " << nam[j] << " =                   " << e
           << " +- " << r;
    G4cout << G4endl;
  }
  G4cout << std::setprecision(4) << "Beam Energy                  " << beamEnergy/GeV 
	 << " GeV" << G4endl;
  if(n_lowe > 0)          G4cout << "Number of events E/E0<0.8    " << n_lowe << G4endl; 
  G4cout<<"=================================================================="<<G4endl;
  G4cout<<G4endl;

  // normalise histograms
  for(G4int i=0; i<nHisto; i++) {
    histo->scale(i,x);
  }

  histo->save();
  if(0 < runID) { return; }

  // Acceptance
  EmAcceptance acc;
  G4bool isStarted = false;
  for (j=0; j<nmax; j++) {

    G4double ltrue = limittrue[j];
    if (ltrue < DBL_MAX) {
      if (!isStarted) {
        acc.BeginOfAcceptance("Crystal Calorimeter",n_evt);
        isStarted = true;
      }
      G4double etrue = edeptrue[j];
      G4double rtrue = rmstrue[j];
      acc.EmAcceptanceGauss("Edep"+nam[j],n_evt,edeptr[j],etrue,rtrue,ltrue);
      acc.EmAcceptanceGauss("Erms"+nam[j],n_evt,ermstr[j],rtrue,rtrue,2.0*ltrue);
    }
  }
  if(isStarted) acc.EndOfAcceptance();

  // atom frequency
  G4cout << "   Z  bremsstrahlung photoeffect  compton    conversion" << G4endl;
  for(j=1; j<93; ++j) {
    G4int n1 = G4int(brem[j]*x);
    G4int n2 = G4int(phot[j]*x);
    G4int n3 = G4int(comp[j]*x);
    G4int n4 = G4int(conv[j]*x);
    if(n1 + n2 + n3 + n4 > 0) {
      G4cout << std::setw(4) << j << std::setw(12) << n1 << std::setw(12) << n2
	     << std::setw(12) << n3 << std::setw(12) << n4 << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfEvent()
{
  n_evt++;

  Eabs1  = 0.0;
  Eabs2  = 0.0;
  Eabs3  = 0.0;
  Eabs4  = 0.0;
  Evertex.clear();
  Nvertex.clear();
  for (G4int i=0; i<25; i++) {
    E[i] = 0.0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfEvent()
{
  G4double e9 = 0.0;
  G4double e25= 0.0;
  for (G4int i=0; i<25; i++) {
    E[i] /= beamEnergy;
    e25 += E[i];
    if( ( 6<=i &&  8>=i) || (11<=i && 13>=i) || (16<=i && 18>=i)) e9 += E[i];
  }

  if(1 < verbose && e25 < 0.8) {
    n_lowe++;
    G4cout << "### in the event# " << n_evt << "  E25= " << e25 << G4endl;
  }

  // compute ratios
  G4double e0 = E[12];
  G4double e19  = 0.0;
  G4double e125 = 0.0;
  G4double e925 = 0.0;
  if(e9 > 0.0) {
    e19 = e0/e9;
    e125 = e0/e25;
    e925 = e9/e25;
    edep[3] += e19;
    erms[3] += e19*e19;
    edep[4] += e125;
    erms[4] += e125*e125;
    edep[5] += e925;
    erms[5] += e925*e925;
    stat[3] += 1;
    stat[4] += 1;
    stat[5] += 1;
  }

  // fill histo
  histo->fill(0,e0,1.0);
  histo->fill(1,e9,1.0);
  histo->fill(2,e25,1.0);
  histo->fill(5,Eabs1,1.0);
  histo->fill(6,Eabs2,1.0);
  histo->fill(7,Eabs3,1.0);
  histo->fill(8,Eabs4,1.0);
  histo->fill(9,G4double(Nvertex.size()),1.0);
  histo->fill(10,e19,1.0);
  histo->fill(11,e125,1.0);
  histo->fill(12,e925,1.0);

  // compute sums
  edep[0] += e0;
  erms[0] += e0*e0;
  edep[1] += e9;
  erms[1] += e9*e9;
  edep[2] += e25;
  erms[2] += e25*e25;

  // trancated mean
  if(limittrue[0] == DBL_MAX || std::abs(e0-edeptrue[0])<rmstrue[0]*limittrue[0]) {
    stat[0] += 1;
    edeptr[0] += e0;
    ermstr[0] += e0*e0;
  }
  if(limittrue[1] == DBL_MAX || std::abs(e9-edeptrue[1])<rmstrue[1]*limittrue[1]) {
    stat[1] += 1;
    edeptr[1] += e9;
    ermstr[1] += e9*e9;
  }
  if(limittrue[2] == DBL_MAX || std::abs(e25-edeptrue[2])<rmstrue[2]*limittrue[2]) {
    stat[2] += 1;
    edeptr[2] += e25;
    ermstr[2] += e25*e25;
  }
  if(nTuple) histo->addRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::ScoreNewTrack(const G4Track* aTrack)
{
  //Save primary parameters
  ResetTrackLength();
  const G4ParticleDefinition* particle = aTrack->GetDefinition();
  const G4DynamicParticle* dynParticle = aTrack->GetDynamicParticle();
  G4String name = particle->GetParticleName();
  G4int pid = aTrack->GetParentID();
  G4double kinE = dynParticle->GetKineticEnergy();
  G4ThreeVector pos = aTrack->GetVertexPosition();

  // primary
  if(0 == pid) {

    beamEnergy = kinE;
    histo->fillTuple("TKIN", kinE/MeV);

    G4double mass = 0.0;
    if(particle) {
      mass = particle->GetPDGMass();
      histo->fillTuple("MASS", mass/MeV);
      histo->fillTuple("CHAR",(particle->GetPDGCharge())/eplus);
      G4double beta = 1.;
	if(mass > 0.) {
          G4double gamma = kinE/mass + 1.;
          beta = std::sqrt(1. - 1./(gamma*gamma));
	}
    }

    G4ThreeVector dir = dynParticle->GetMomentumDirection();
    if(1 < verbose) {
      G4cout << "TrackingAction: Primary kinE(MeV)= " << kinE/MeV
           << "; m(MeV)= " << mass/MeV
           << "; pos= " << pos << ";  dir= " << dir << G4endl;
    }

    // secondary
  } else {
    const G4VProcess* proc = aTrack->GetCreatorProcess();
    G4int type = proc->GetProcessSubType();
    if(type == fBremsstrahlung) {
      const G4Element* elm = static_cast<const G4VEnergyLossProcess*>(proc)->GetCurrentElement();
      if(elm) {
	G4int Z = G4int(elm->GetZ());
        if(Z > 0 && Z < 93) { brem[Z] += 1.0; }
      }
    } else if(type == fPhotoElectricEffect) {
      const G4Element* elm = static_cast<const G4VEmProcess*>(proc)->GetCurrentElement();
      if(elm) {
	G4int Z = G4int(elm->GetZ());
        if(Z > 0 && Z < 93) { phot[Z] += 1.0; }
      }
    } else if(type == fGammaConversion) {
      const G4Element* elm = static_cast<const G4VEmProcess*>(proc)->GetCurrentElement();
      if(elm) {
	G4int Z = G4int(elm->GetZ());
        if(Z > 0 && Z < 93) { conv[Z] += 1.0; }
      }
    } else if(type == fComptonScattering) {
      const G4Element* elm = static_cast<const G4VEmProcess*>(proc)->GetCurrentElement();
      if(elm) {
	G4int Z = G4int(elm->GetZ());
        if(Z > 0 && Z < 93) { comp[Z] += 1.0; }
      }
    }

    // delta-electron
    if (0 < pid && "e-" == name) {
      if(1 < verbose) {
	G4cout << "TrackingAction: Secondary electron " << G4endl;
      }
      AddDeltaElectron(dynParticle);

    } else if (0 < pid && "e+" == name) {
      if(1 < verbose) {
	G4cout << "TrackingAction: Secondary positron " << G4endl;
      }
      AddPositron(dynParticle);

    } else if (0 < pid && "gamma" == name) {
      if(1 < verbose) {
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
  if(1 < verbose) {
    G4cout << "HistoManager::AddEnergy: e(keV)= " << edep/keV
           << "; volIdx= " << volIndex
           << "; copyNo= " << copyNo
           << G4endl;
  }
  if(0 == volIndex) {
    E[copyNo] += edep;
  } else if (1 == volIndex) {
    Eabs1 += edep;
  } else if (2 == volIndex) {
    Eabs2 += edep;
  } else if (3 == volIndex) {
    Eabs3 += edep;
  } else if (4 == volIndex) {
    Eabs4 += edep;
  } else if (5 == volIndex) {
    G4int n = Nvertex.size();
    G4bool newPad = true;
    if (n > 0) {
      for(G4int i=0; i<n; i++) {
        if (copyNo == Nvertex[i]) {
          newPad = false;
          Evertex[i] += edep;
          break;
        }
      }
    }
    if(newPad) {
      Nvertex.push_back(copyNo);
      Evertex.push_back(edep);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddDeltaElectron(const G4DynamicParticle* elec)
{
  G4double e = elec->GetKineticEnergy()/MeV;
  if(e > 0.0) n_elec++;
  histo->fill(3,e,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddPhoton(const G4DynamicParticle* ph)
{
  G4double e = ph->GetKineticEnergy()/MeV;
  if(e > 0.0) n_gam++;
  histo->fill(4,e,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetEdepAndRMS(G4int i, G4ThreeVector val)
{
  if(i<nmax && i>=0) {
    if(val[0] > 0.0) edeptrue[i] = val[0];
    if(val[1] > 0.0) rmstrue[i] = val[1];
    if(val[2] > 0.0) limittrue[i] = val[2];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

