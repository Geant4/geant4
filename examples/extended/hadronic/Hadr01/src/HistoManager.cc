//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "Histo.hh"
#include "G4NistManager.hh"
#include "globals.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"

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
  verbose=  0;
  nEvt1  = -1;
  nEvt2  = -1;

  nTuple = true;

  maxNeutronEnergy  = 20.0*MeV;
  maxHadronEnergy   = 100.0*MeV;
  maxElectronEnergy = 20.0*MeV;
  maxGammaEnergy    = 20.0*MeV;
  nBinsE    = 100;
  nBinXY    = 20;
  n_XY      = 20;
  nHisto    = 17;
  absWidth  = 1.*mm;
  gapWidth  = 0.*mm;
  absRadius = 10.*cm;
  n_abs     = 300;

  Z0 = 8.0;
  A0 = 14.0;
  histo = new Histo(verbose);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete histo;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::bookHisto()
{
  absLength = absWidth*(G4double)n_abs;
  absZ0     = -0.5*absLength;

  histo->add1D("10",
    "Energy deposition (MeV/mm/event) in phantom",n_abs,0.0,absLength,MeV/mm);

  histo->add1D("11",
    "Energy (MeV) of neutrons",nBinsE,0.0,maxHadronEnergy,MeV);

  histo->add1D("12",
    "Energy (MeV) of delta-electrons",nBinsE,0.0,maxElectronEnergy,MeV);

  histo->add1D("13",
    "Energy (MeV) of gammas",nBinsE,0.0,maxGammaEnergy,MeV);

  histo->add1D("14",
    "Energy of side-leaked neutrons (MeV)",nBinsE,0.0,maxNeutronEnergy,MeV);

  histo->add1D("15",
    "Energy of forward-leaked neutrons (MeV)",nBinsE,0.0,maxNeutronEnergy,MeV);

  histo->add1D("16",
    "Energy of backward-leaked neutrons (MeV)",nBinsE,0.0,maxNeutronEnergy,MeV);

  histo->add1D("17",
    "Energy deposition over X in the Bragg Peak (MeV/mm/event) in phantom",
	      n_XY,-absRadius,absRadius,MeV/mm);

  histo->add1D("18",
    "Energy deposition over Y in the Bragg Peak (MeV/mm/event) in phantom",
	      n_XY,-absRadius,absRadius,MeV/mm);

  histo->add1D("19",
    "Isotope production (1/cm/event)",n_abs/10,0.0,absLength,cm);

  histo->add1D("20",
    "Energy of protons (MeV)",nBinsE,0.0,maxHadronEnergy,MeV);

  histo->add1D("21",
    "Energy of pions (MeV)",nBinsE,0.0,maxHadronEnergy,MeV);

  histo->add1D("22",
    "Log10 Energy (MeV) of gammas",100,-5.,5.,1.0);

  histo->add1D("23",
    "Log10 Energy (MeV) of charged pions",100,-4.,6.,1.0);

  histo->add1D("24",
    "Log10 Energy (MeV) of protons",100,-5.,5.,1.0);

  histo->add1D("25",
    "Log10 Energy (MeV) of neutrons",100,-5.,5.,1.0);

  histo->add1D("26",
    "Log10 Energy (MeV) of pi0",100,-4.,6.,1.0);
  /*
 
  histo->add1D("31",
    "Isotope Mass Difference (MeV) NIST - G4IonTable",100,-0.5,99.5,1.0);
  histo->add1D("32",
    "Atomic Mass Difference (MeV) GHAD - G4IonTable",100,-0.5,99.5,1.0);
  histo->add1D("33",
    "Atomic Mass Difference (MeV) GHAD - G4IonTable - ZxMe",100,-0.5,99.5,1.0);
  */
  if(nTuple) {
    histo->addTuple( "100", "tuple.paw","float Z, A, nist, ghad" );

    // histo->addTuple( "100", "tuple.paw","float zz, Z, A" );

    G4cout<<" addTuple (HistoManager). "<<G4endl;
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
  n_prot_leak = 0;
  n_aprot_leak= 0;
  n_pos_hadr = 0;
  n_neg_hadr = 0;
  n_ions     = 0;
  n_neutron  = 0;
  n_neu_forw = 0;
  n_neu_leak = 0;
  n_neu_back = 0;
  endX2  = 0.0;
  endY2  = 0.0;
  endZ2  = 0.0;
  endZ   = 0.0;
  trackLength = 0.0;

  bookHisto();
  histo->book();

  if(verbose > 0) 
  {
    G4cout << "HistoManager: Histograms are booked and run has been started"
           <<G4endl<<"  BeginOfRun (After histo->book)"<< G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun()
{

  G4cout << "Histo: End of run actions are started" << G4endl;

  // average

  G4cout<<"========================================================"<<G4endl;
  G4double x = (G4double)n_evt;
  if(n_evt > 0) x = 1.0/x;
  G4double xe = x*(G4double)n_elec;
  G4double xg = x*(G4double)n_gam;
  G4double xp = x*(G4double)n_posit;
  G4double xs = x*(G4double)n_step;
  G4double xn = x*(G4double)n_neutron;
  G4double xnf = x*(G4double)n_neu_forw;
  G4double xnb = x*(G4double)n_neu_leak;
  G4double xnbw= x*(G4double)n_neu_back;
  G4double xpl = x*(G4double)n_prot_leak;
  G4double xal = x*(G4double)n_aprot_leak;
  G4double xph = x*(G4double)n_pos_hadr;
  G4double xnh = x*(G4double)n_neg_hadr;
  G4double xio = x*(G4double)n_ions;
  endZ  *= x;
  endZ2 *= x;
  G4cout                         << "Beam particle                        "
				 << primaryDef->GetParticleName() <<G4endl;
  G4cout                         << "Beam Energy(MeV)                     " 
				 << primaryKineticEnergy/MeV <<G4endl;
  G4cout                         << "Number of events                     " << n_evt <<G4endl;
  G4cout << std::setprecision(4) << "Average number of e-                 " << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma              " << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+                 " << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of neutrons           " << xn << G4endl;
  G4cout << std::setprecision(4) << "Average number of forward neutrons   " << xnf << G4endl;
  G4cout << std::setprecision(4) << "Average number of back neutrons      " << xnb << G4endl;
  G4cout << std::setprecision(4) << "Average number of backword neutrons  " << xnbw << G4endl;
  G4cout << std::setprecision(4) << "Average number of proton leak        " << xpl << G4endl;
  G4cout << std::setprecision(4) << "Average number of anti proton leak   " << xal << G4endl;
  G4cout << std::setprecision(4) << "Average number of positive hadrons   " << xph << G4endl;
  G4cout << std::setprecision(4) << "Average number of negative hadrons   " << xnh << G4endl;
  G4cout << std::setprecision(4) << "Average number of ions               " << xio << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps              " << xs << G4endl;
  G4cout << std::setprecision(4) << "Average track length (mm)            " << x*trackLength/mm << G4endl;
  G4cout << std::setprecision(4) << "End point: sigma X (mm)              " << std::sqrt(x*endX2) << G4endl;
  G4cout << std::setprecision(4) << "End point: sigma Y (mm)              " << std::sqrt(x*endY2) << G4endl;
  G4cout << std::setprecision(4) << "End point: sigma Z (mm)              " 
	 << std::sqrt(endZ2 - endZ*endZ) << G4endl;
  G4cout << std::setprecision(4) << "End point: position Z (mm)           " << endZ << G4endl;
  G4cout<<"========================================================"<<G4endl;
  G4cout<<G4endl;

  // normalise histograms
  for(G4int i=0; i<nHisto; i++) {
    histo->scale(i,x);
  }

  histo->print(0);
  IsotopeStudy();
  histo->save();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfEvent()
{
  n_evt++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfEvent()
{
  if(nTuple) histo->addRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::ScoreNewTrack(const G4Track* aTrack)
{
  //Save primary parameters
  const G4ParticleDefinition* particle = aTrack->GetDefinition();
  const G4DynamicParticle* dynParticle = aTrack->GetDynamicParticle();
  G4String name = particle->GetParticleName();
  G4int pid = aTrack->GetParentID();
  G4double kinE = dynParticle->GetKineticEnergy();
  G4ThreeVector pos = aTrack->GetVertexPosition();

  if(0 == pid) {

    primaryKineticEnergy = kinE;
    primaryDef = particle;
    G4ThreeVector dir = dynParticle->GetMomentumDirection();
    if(1 < verbose) {
      G4cout << "TrackingAction: Primary kinE(MeV)= " << kinE/MeV
           << "; m(MeV)= " << dynParticle->GetMass()/MeV
           << "; pos= " << pos << ";  dir= " << dir << G4endl;
    }

    // delta-electron
  } else {
    if ("e-" == name) {
      if(1 < verbose) {
        G4cout << "TrackingAction: Secondary electron "
               << " E= " << kinE << G4endl;
     }
      AddElectron(dynParticle);

    } else if ("e+" == name) {
      if(1 < verbose) {
        G4cout << "TrackingAction: Secondary positron "
               << " E= " << kinE << G4endl;
      }
      AddPositron(dynParticle);

    } else if ("gamma" == name) {
      if(1 < verbose) {
        G4cout << "TrackingAction: Secondary gamma; parentID= " << pid
               << " E= " << kinE << G4endl;
      }
      AddPhoton(dynParticle);

    } else if ("neutron" == name) {
      if(1 < verbose) {
        G4cout << "TrackingAction: Secondary neutron; parentID= " << pid
               << " E= " << kinE << G4endl;
      }
      AddNeutron(dynParticle);

    } else if ("proton" == name) {
      if(1 < verbose) {
        G4cout << "TrackingAction: Secondary neutron; parentID= " << pid
               << " E= " << kinE << G4endl;
      }
      AddProton(dynParticle);

    } else {
      if (particle->GetParticleType() == "nucleus") {
	if(1 < verbose) {
	  G4cout << "TrackingAction: Secondary ion; parentID= " << pid
		 << " E= " << kinE << G4endl;
	}
	AddIon(dynParticle, pos);
      } else {
	AddHadron(dynParticle);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddEnergy(G4double edep, G4double step, const G4ThreeVector& v)
{
  G4double x = v.x();
  G4double y = v.y();
  G4double z = v.z() - absZ0;
  if(1 < verbose) {
    G4cout << "HistoManager::AddEnergy: e(keV)= " << edep/keV
           << "; step(mm)= " << step/mm
           << "; z(mm)= " << z/mm
           << G4endl;
  }
  histo->fill(0,z,edep);
  G4int bin = G4int(z/absWidth);
  if(bin == nBinXY) {
    histo->fill(7,x,edep);
    histo->fill(8,y,edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddNeutron(const G4DynamicParticle* neu)
{
  G4double e = neu->GetKineticEnergy();
  G4double mom = neu->GetTotalMomentum();
  if(e > 0.0) n_neutron++;
  //  histo->fill(1,e,1.0);
  histo->fill(15,std::log10(e/MeV),1.0);
  if(mom > 0.0)
    histo->fill(1,std::log10(e/MeV),1.0/mom);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddProton(const G4DynamicParticle* p)
{
  G4double e = p->GetKineticEnergy();
  G4double mom = p->GetTotalMomentum();
  if(e > 0.0) n_pos_hadr++;
  histo->fill(14,std::log10(e/MeV),1.0);
  if(mom > 0.0)
    histo->fill(10,std::log10(e/MeV),1.0/mom);
}        

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddPhoton(const G4DynamicParticle* p)
{
  G4double e = p->GetKineticEnergy();
  n_gam++;
  histo->fill(3,e,1.0);
  histo->fill(12,std::log10(e/MeV),1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddHadron(const G4DynamicParticle* p)
{
  G4double e = p->GetKineticEnergy();
  G4double m = p->GetMass()/MeV;
  G4double mom = p->GetTotalMomentum();
  if(m < 150. && m > 138. ) {
    histo->fill(13,std::log10(e/MeV),1.0);
    if(mom > 0.0)
      histo->fill(11,std::log10(e/MeV),1.0/mom);
  } else if(m < 138. && m > 130. ) {
    histo->fill(16,std::log10(e/MeV),1.0);
    if(mom > 0.0)
      histo->fill(11,std::log10(e/MeV),1.0/mom);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddElectron(const G4DynamicParticle* elec)
{
  G4double e = elec->GetKineticEnergy();
  n_elec++;
  histo->fill(2,e,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddPositron(const G4DynamicParticle*)
{
  n_posit++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddLeakingNeutron(G4double e, const G4ThreeVector& pos)
{
  G4double x = pos.x();
  G4double y = pos.y();
  G4double z = pos.z();
  G4double r = std::sqrt(x*x + y*y);
  if(r > absRadius) {
    n_neu_leak++;
    histo->fill(4,e,1.0);
  } else if (z > 0.0) {
    n_neu_forw++;
    histo->fill(5,e,1.0);
  } else {
    n_neu_back++;
    histo->fill(6,e,1.0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddLeakingHadron(G4double e, const G4ThreeVector& pos, 
				    const G4ParticleDefinition* p)
{
  G4double x = pos.x();
  G4double y = pos.y();
  //  G4double z = pos.z();
  G4double r = std::sqrt(x*x + y*y);
  G4double q = p->GetPDGCharge();
  if(r > absRadius) {
    if(q > 0.0) {
      histo->fill(10,e,1.0);
      if(p == G4Proton::Proton()) n_prot_leak++;
    } else if (q < 0.0) {
      histo->fill(11,e,1.0);
      if(p == G4AntiProton::AntiProton()) n_aprot_leak++;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::AddIon(const G4DynamicParticle* ph, const G4ThreeVector& pos)
{
  n_ions++;
  G4double z = pos.z() - absZ0;
  //  histo->fillTuple("zz", z/mm);
  const G4ParticleDefinition* pd = ph->GetDefinition();
  G4double Z = pd->GetPDGCharge()/eplus;
  G4double A = pd->GetBaryonNumber();
  // histo->fillTuple("Z", Z);
  // histo->fillTuple("A", A);
  // histo->addRow();
  if(std::fabs(Z-Z0)<0.5 && std::abs(A-A0)<0.5) histo->fill(9,z,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetEndPoint(const G4ThreeVector& pos)
{
  G4double x = pos.x();
  G4double y = pos.y();
  G4double z = pos.z() - absZ0;
  endX2 += x*x;
  endY2 += y*y;
  endZ  += z;
  endZ2 += z*z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetVerbose(G4int val)        
{
  verbose = val; 
  histo->setVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::IsotopeStudy()
{
  G4cout <<"=========================================================="<<G4endl;
  G4NistManager* man = G4NistManager::Instance();
  G4IonTable itable;
  G4Nucleus  ghad;
  G4int i,j;
  G4double p[5] = {3.0, 5.0, 8.0, 12.0, 15.};
  G4double m2   = man->GetIsotopeMass(73,181)*amu_c2;
  G4double m1   = proton_mass_c2;
  G4double st   = 0.5;
  G4double ct   = sqrt(1.0 - st*st);
  for(i=0; i<5; i++) {
    G4double p1 = p[i]*GeV;
    G4double e1 = sqrt(p1*p1 + m1*m1);
    G4double p2 = 2.0*p1*(e1 + m2)*m2*ct/(p1*p1*st*st + m1*m1 + m2*m2 + 2.0*e1*m2);
    G4double e2 = sqrt(p2*p2 + m2*m2);
    G4double e3 = e1 + m2 - e2; 
    G4double p3x= p1 - p2*ct; 
    G4double p3y= p2*st;
    G4double m3 = sqrt(e3*e3 - p3x*p3x - p3y*p3y) - m1; 
    G4cout << " p(GeV)= " << p1/GeV << "   At 45deg p2(MeV)= "<< p2/MeV 
	   << "  e2(MeV)= " << e2 - m2 << "  del= " << m3 <<G4endl; 
  }
  G4cout <<"=========================================================="<<G4endl;
  for(i=1; i<93; i++) {
    G4double Z = G4double(i);
    G4double amax = 0.0;
    G4double nist0= 0.0;
    G4double nist = 0.0;
    G4double gh   = 0.0;
    for(j=1; j<270; j++) {
      G4double mNist = man->GetIsotopeMass(i,j)*amu_c2;
      G4double A = G4double(j);
      if(mNist > 0.0) {
	G4double aNist = man->GetIsotopeAbundance(i,j)*100.;
        G4double mPart = itable.GetNucleusMass(i,j);
        G4double mGhad = ghad.AtomicMass(A, Z) -  mPart;
        nist0 = mNist - mPart;
 
	G4cout << "Z= " << i << " A= " << j << " aband(%)= " << aNist
	       << " errP(MeV)= " << nist0 << " errG(MeV)= " << mGhad
	       << G4endl; 
	histo->fillTuple("Z", Z);
	histo->fillTuple("A", A);
	histo->fillTuple("nist", nist0);
	histo->fillTuple("ghad", mGhad);
	histo->addRow();
	if(aNist > amax) {
          amax = aNist;
          nist = nist0;
          gh   = mGhad;
	}
      }
    }
    histo->fill(17,Z,nist);
    histo->fill(18,Z,gh); 
    histo->fill(19,Z,gh-Z*electron_mass_c2); 
  }
  G4cout <<"=========================================================="<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

