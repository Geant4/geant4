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
  nHisto = 10;

  for(G4int i=0; i<nmax; i++) {
    edep.push_back(0.0);
    erms.push_back(0.0);
    edeptrue.push_back(1.0);
    rmstrue.push_back(1.0);
    limittrue.push_back(DBL_MAX);
  }

  histo->add1D("10",
    "Energy deposit (MeV) in central crystal",nBinsED,0.0,beamEnergy,MeV);

  histo->add1D("11",
    "Energy deposit (MeV) in 3x3",nBinsED,0.0,beamEnergy,MeV);

  histo->add1D("12",
    "Energy deposit (MeV) in 5x5",nBinsED,0.0,beamEnergy,MeV);

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

  if(nTuple) {
    histo->addTuple( "100", "Dose deposite","float r, z, e" );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  n_evt  = 0;
  n_elec = 0;
  n_posit= 0;
  n_gam  = 0;
  n_step = 0;
  histo->book();

  if(verbose > 0) {
    G4cout << "HistoManager: Histograms are booked and run has been started"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun()
{

  G4cout << "HistoManager: End of run actions are started" << G4endl;
  G4String nam[3] = {"1x1", "3x3", "5x5"};

  // average

  G4cout<<"========================================================"<<G4endl;
  G4double x = (G4double)n_evt;
  if(n_evt > 0) x = 1.0/x;
  G4int j;
  for(j=0; j<nmax; j++) {
    edep[j] *= x/beamEnergy;
    G4double y = erms[j]*x/(beamEnergy*beamEnergy) - edep[j]*edep[j];
    if(y < 0.0) y = 0.0;
    erms[j] = std::sqrt(y);
  }
  G4double xe = x*(G4double)n_elec;
  G4double xg = x*(G4double)n_gam;
  G4double xp = x*(G4double)n_posit;
  G4double xs = x*(G4double)n_step;

  G4double f = 100.*std::sqrt(beamEnergy/GeV);

  G4cout                         << "Number of events             " << n_evt <<G4endl;
  G4cout << std::setprecision(4) << "Average number of e-         " << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma      " << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+         " << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps      " << xs << G4endl;
  for(j=0; j<3; j++) {
    G4cout << std::setprecision(4) << "Edep " << nam[j] << " =                   " << edep[j]
           << " +- " << erms[j]*std::sqrt(x);
    if(edep[j] > 0.0)
      G4cout << "  res=  " << f*erms[j]/edep[j] << " %";
    G4cout << G4endl;
  }
  G4cout<<"========================================================"<<G4endl;
  G4cout<<G4endl;

  // normalise histograms
  for(G4int i=0; i<nHisto; i++) {
    histo->scale(i,x);
  }

  histo->save();

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
      acc.EmAcceptanceGauss("Edep"+nam[j],n_evt,edep[j],etrue,rtrue,ltrue);
      acc.EmAcceptanceGauss("Erms"+nam[j],n_evt,erms[j],rtrue,rtrue,2.0*ltrue);
    }
  }
  if(isStarted) acc.EndOfAcceptance();
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
  for (int i=0; i<25; i++) {
    E[i] = 0.0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfEvent()
{
  G4double e9 = 0.0;
  G4double e25= 0.0;
  for (int i=0; i<25; i++) {
    e25 += E[i];
    if( ( 6<=i &&  8>=i) || (11<=i && 13>=i) || (16<=i && 18>=i)) e9 += E[i];
  }
  histo->fill(0,E[12],1.0);
  histo->fill(1,e9,1.0);
  histo->fill(2,e25,1.0);
  histo->fill(5,Eabs1,1.0);
  histo->fill(6,Eabs2,1.0);
  histo->fill(7,Eabs3,1.0);
  histo->fill(8,Eabs4,1.0);
  float nn = (double)(Nvertex.size());
  histo->fill(9,nn,1.0);

  edep[0] += E[12];
  erms[0] += E[12]*E[12];
  edep[1] += e9;
  erms[1] += e9*e9;
  edep[2] += e25;
  erms[2] += e25*e25;

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

  if(0 == pid) {

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

    // delta-electron
  } else if (0 < pid && "e-" == name) {
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

