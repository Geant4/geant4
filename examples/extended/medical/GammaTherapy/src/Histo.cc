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
// $Id: Histo.cc,v 1.10 2010-10-26 12:05:14 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4Neutron.hh"
#include <iomanip>

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

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
  verbose   = 1;
  histName  = G4String("histo");
  histType  = G4String("root");
  nHisto    = 10;
  nHisto1   = 10;
  maxEnergy = 50.0*MeV;
  nTuple = false;
  nBinsZ = 60;
  nBinsR = 80;
  nBinsE = 200;
  absorberZ = 300.*mm;
  absorberR = 200.*mm;
  scoreZ    = 100.*mm;

  gamma    = G4Gamma::Gamma();
  electron = G4Electron::Electron();
  positron = G4Positron::Positron();
  neutron  = G4Neutron::Neutron();

#ifdef G4ANALYSIS_USE
  af   = 0;
  tree = 0;
  ntup = 0;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::~Histo()
{
#ifdef G4ANALYSIS_USE
  delete af;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::BeginOfHisto()
{
  G4cout << "### Histo start initialisation nHisto= " << nHisto << G4endl;

  n_evt  = 0;
  n_elec = 0;
  n_posit= 0;
  n_gam  = 0;
  n_step = 0;
  n_gam_ph= 0;
  n_gam_tar= 0;
  n_e_tar = 0;
  n_e_ph  = 0;
  n_step_target = 0;
  n_neutron = 0;
  sumR = 0.0;
  if(nBinsR>1000) SetNumberDivR(40);

  stepZ = absorberZ/(G4double)nBinsZ;
  stepR = absorberR/(G4double)nBinsR;
  stepE = maxEnergy/(G4double)nBinsE;
  nScoreBin = (G4int)(scoreZ/stepZ + 0.5);

  G4cout << "   "<< nBinsR << " bins R   stepR= " << stepR/mm << " mm " << G4endl;
  G4cout << "   "<< nBinsZ << " bins Z   stepZ= " << stepZ/mm << " mm " << G4endl;
  G4cout << "   "<< nBinsE << " bins E   stepE= " << stepE/MeV << " MeV " << G4endl;
  G4cout << "   "<< nScoreBin << "th bin in Z is used for R distribution" << G4endl;

  G4int i;
  G4double r1 = 0.0;
  G4double r2 = stepR;
  volumeR.clear();
  for(i=0; i<nBinsR; i++) {
    volumeR.push_back(cm*cm/(pi*(r2*r2 - r1*r1)));
    r1 = r2;
    r2 += stepR;
  }
  for(i=0; i<nBinsE; i++) {
    gammaE.push_back(0.0);
  }

  bookHisto();

  if(verbose > 0) {
    G4cout << "Histo: Histograms are booked and run has been started"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::EndOfHisto()
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
  G4double xph= x*(G4double)n_gam_ph;
  G4double xes= x*(G4double)n_step_target;
  G4double xgt= x*(G4double)n_gam_tar;
  G4double xet= x*(G4double)n_e_tar;
  G4double xphe= x*(G4double)n_e_ph;
  G4double xne= x*(G4double)n_neutron;
  G4cout                    << "Number of events                                  " << n_evt <<G4endl;
  G4cout << std::setprecision(4) << "Average number of e-                         " << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma                      " << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+                         " << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of neutrons                   " << xne << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps in absorber          " << xs << G4endl;
  G4cout << std::setprecision(4) << "Average number of e- steps in target         " << xes << G4endl;
  G4cout << std::setprecision(4) << "Average number of g  produced in the target  " << xgt << G4endl;
  G4cout << std::setprecision(4) << "Average number of e- produced in the target  " << xet << G4endl;
  G4cout << std::setprecision(4) << "Average number of g produced in the phantom  " << xph << G4endl;
  G4cout << std::setprecision(4) << "Average number of e- produced in the phantom " << xphe << G4endl;
  G4cout << std::setprecision(4) << "Total gamma fluence in front of phantom      " << x*sumR/MeV
         << " MeV " << G4endl;
  G4cout<<"========================================================"<<G4endl;
  G4cout<<G4endl;
  G4cout<<G4endl;

#ifdef G4ANALYSIS_USE
  if(tree) {
    // normalise histograms
    for(G4int i=0; i<nHisto1; i++) {
      histo[i]->scale(x);
    }
    G4double nr = histo[0]->binHeight(0);
    if(nr > 0.0) nr = 1./nr;
    histo[0]->scale(nr);

    nr = (histo[1]->sumAllBinHeights())*stepR;
    if(nr > 0.0) nr = 1./nr;
    histo[1]->scale(nr);

    histo[3]->scale(1000.0*cm3/(pi*absorberR*absorberR*stepZ));
    histo[4]->scale(1000.0*cm3*volumeR[0]/stepZ);

    // Write histogram file
    if(0 < nHisto) {
      tree->commit();
      G4cout << "Histograms and Ntuples are saved" << G4endl;
    }
    tree->close();
    delete tree;
    tree = 0;
    G4cout << "Tree is closed" << G4endl;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::SaveEvent()
{
#ifdef G4ANALYSIS_USE
  if(ntup) {
    ntup->addRow();
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::SaveToTuple(const G4String& parname, G4double val)
{
  if(2 < verbose) G4cout << "Save to tuple " << parname << "   " << val << G4endl;
#ifdef G4ANALYSIS_USE
  if(ntup) ntup->fill( ntup->findColumn(parname), (float)val);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::SaveToTuple(const G4String& parname,G4double val, G4double)
{
  if(2 < verbose) G4cout << "Save to tuple " << parname << "   " << val << G4endl;
#ifdef G4ANALYSIS_USE
  if(ntup) ntup->fill( ntup->findColumn(parname), (float)val);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::bookHisto()
{

#ifdef G4ANALYSIS_USE
  G4cout << "### Histo books " << nHisto << " histograms " << G4endl;
  // Creating the analysis factory
  if(!af) af = AIDA_createAnalysisFactory();

  // Creating the tree factory
  AIDA::ITreeFactory* tf = af->createTreeFactory();

  // Creating a tree mapped to a new hbook file.
  G4String tt = histType;
  G4String nn = histName + "." + histType;
  if(histType == "xml" || histType == "XML" || histType == "aida" || 
     histType == "AIDA") {
    tt = "xml";
    nn = histName + ".aida";
  }

  tree = tf->create(nn,tt,false,true, "");
  if(tree) {
    G4cout << "Tree store : " << tree->storeName() << G4endl;
  } else {
    G4cout << "Fail to open tree store " << nn << G4endl;
    return;
  }
  delete tf;
  histo.resize(nHisto1);

  // Creating a histogram factory, whose histograms will be handled by the tree
  AIDA::IHistogramFactory* hf = af->createHistogramFactory( *tree );

  // Creating an 1-dimensional histograms in the root directory of the tree

  histo[0] = hf->createHistogram1D("10",
    "Energy deposit at radius (mm) normalised on 1st channel",nBinsR,0.,absorberR/mm);

  histo[1] = hf->createHistogram1D("11",
    "Energy deposit at radius (mm) normalised to integral",nBinsR,0.,absorberR/mm);

  histo[2] = hf->createHistogram1D("12",
    "Energy deposit (MeV/kg/electron) at radius (mm)",nBinsR,0.,absorberR/mm);

  histo[3] = hf->createHistogram1D("13",
    "Energy profile (MeV/kg/electron) over Z (mm)",nBinsZ,0.,absorberZ/mm);

  histo[4] = hf->createHistogram1D("14",
    "Energy profile (MeV/kg/electron) over Z (mm) at Central Voxel",nBinsZ,0.,absorberZ/mm);

  histo[5] = hf->createHistogram1D("15",
    "Energy (MeV) of gamma produced in the target",nBinsE,0.,maxEnergy/MeV);

  histo[6] = hf->createHistogram1D("16",
    "Energy (MeV) of gamma before phantom",nBinsE,0.,maxEnergy/MeV);

  histo[7] = hf->createHistogram1D("17",
    "Energy (MeV) of electrons produced in phantom",nBinsE,0.,maxEnergy/MeV);

  histo[8] = hf->createHistogram1D("18",
    "Energy (MeV) of electrons produced in target",nBinsE,0.,maxEnergy/MeV);

  histo[9] = hf->createHistogram1D("19",
    "Gamma Energy Fluence (MeV/cm2) at radius(mm) in front of phantom",nBinsR,0.,absorberR/mm);

  // Creating a tuple factory, whose tuples will be handled by the tree
  AIDA::ITupleFactory* tpf = af->createTupleFactory( *tree );

  // If using Anaphe HBOOK implementation, there is a limitation on the
  // length of the variable names in a ntuple
  if(nTuple) {
     ntup = tpf->create( "100", "Dose deposite","float r, z, e" );
  }
  delete hf;
  delete tpf;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddDeltaElectron(const G4DynamicParticle* elec)
{
  G4double e = elec->GetKineticEnergy()/MeV;
  if(e > 0.0) { ++n_elec; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddPhoton(const G4DynamicParticle* ph)
{
  G4double e = ph->GetKineticEnergy()/MeV;
  if(e > 0.0) { ++n_gam; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddTargetPhoton(const G4DynamicParticle* ph)
{
  G4double e = ph->GetKineticEnergy()/MeV;
  if(e > 0.0) { ++n_gam_tar; }
#ifdef G4ANALYSIS_USE
  if(tree) histo[5]->fill(e,1.0);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddPhantomPhoton(const G4DynamicParticle* ph)
{
  G4double e = ph->GetKineticEnergy()/MeV;
  if(e > 0.0) { ++n_gam_ph; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddTargetElectron(const G4DynamicParticle* ph)
{
  G4double e = ph->GetKineticEnergy()/MeV;
  if(e > 0.0) { ++n_e_tar; }
#ifdef G4ANALYSIS_USE
  if(tree) histo[8]->fill(e,1.0);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddPhantomElectron(const G4DynamicParticle* ph)
{
  G4double e = ph->GetKineticEnergy()/MeV;
  if(e > 0.0) { ++n_e_ph; }
#ifdef G4ANALYSIS_USE
  if(tree) histo[7]->fill(e,1.0);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::ScoreNewTrack(const G4Track* aTrack)
{
  //Save primary parameters

  ResetTrackLength();
  const G4ParticleDefinition* particle = aTrack->GetParticleDefinition();
  G4int pid = aTrack->GetParentID();
  G4double kinE = aTrack->GetKineticEnergy();
  G4ThreeVector pos = aTrack->GetVertexPosition();
  G4VPhysicalVolume* pv = aTrack->GetVolume();
  const G4DynamicParticle* dp = aTrack->GetDynamicParticle();

  if(0 == pid) {

    SaveToTuple("TKIN", kinE/MeV);

    G4double mass = 0.0;
    if(particle) {
      mass = particle->GetPDGMass();
      SaveToTuple("MASS", mass/MeV);
      SaveToTuple("CHAR",(particle->GetPDGCharge())/eplus);
    }

    G4ThreeVector dir = aTrack->GetMomentumDirection();
    if(1 < verbose) {
      G4cout << "TrackingAction: Primary kinE(MeV)= " << kinE/MeV
           << "; m(MeV)= " << mass/MeV
           << "; pos= " << pos << ";  dir= " << dir << G4endl;
    }

    // delta-electron
  } else if (0 < pid && particle == electron) {
    if(1 < verbose) {
      G4cout << "TrackingAction: Secondary electron " << G4endl;
    }
    AddDeltaElectron(dp);
    if(pv == phantom)                       { AddPhantomElectron(dp); }
    else if(pv == target1 || pv == target2) { AddTargetElectron(dp); }

  } else if (0 < pid && particle == positron) {
    if(1 < verbose) {
      G4cout << "TrackingAction: Secondary positron " << G4endl;
    }
    AddPositron(dp);

  } else if (0 < pid && particle == gamma) {
    if(1 < verbose) {
      G4cout << "TrackingAction: Secondary gamma; parentID= " << pid
             << " E= " << aTrack->GetKineticEnergy() << G4endl;
    }
    AddPhoton(dp);
    if(pv == phantom)                       { AddPhantomPhoton(dp); }
    else if(pv == target1 || pv == target2) { AddTargetPhoton(dp); }

  } else if (0 < pid && particle == neutron) {
    n_neutron++;
    if(1 < verbose) {
      G4cout << "TrackingAction: Secondary neutron; parentID= " << pid
             << " E= " << aTrack->GetKineticEnergy() << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddGamma(G4double e, G4double r)
{
  e /= MeV;
  sumR += e;
  G4int bin = (G4int)(e/stepE);
  if(bin >= nBinsE) { bin = nBinsE-1; }
  gammaE[bin] += e;
  G4int bin1 = (G4int)(r/stepR);
  if(bin1 >= nBinsR) { bin1 = nBinsR-1; }
#ifdef G4ANALYSIS_USE
  if(tree) {
    histo[6]->fill(e,1.0);
    histo[9]->fill(r,e*volumeR[bin1]);
  }
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::AddStep(G4double edep, G4double r1, G4double z1, G4double r2, G4double z2,
                                   G4double r0, G4double z0)
{
  n_step++;
  G4int nzbin = (G4int)(z0/stepZ);
  if(verbose > 1) {
    G4cout << "Histo: edep(MeV)= " << edep/MeV << " at binz= " << nzbin
           << " r1= " << r1 << " z1= " << z1
           << " r2= " << r2 << " z2= " << z2
           << " r0= " << r0 << " z0= " << z0
	   << G4endl;
  }
  if(nzbin == nScoreBin) {
#ifdef G4ANALYSIS_USE
    if(tree) {
      G4int bin  = (G4int)(r0/stepR);
      if(bin >= nBinsR) bin = nBinsR-1;
      double w    = edep*volumeR[bin];
      histo[0]->fill(r0,w);
      histo[1]->fill(r0,w);
      histo[2]->fill(r0,w);
    }
#endif
  }
  G4int bin1 = (G4int)(z1/stepZ);
  if(bin1 >= nBinsZ) bin1 = nBinsZ-1;
  G4int bin2 = (G4int)(z2/stepZ);
  if(bin2 >= nBinsZ) bin2 = nBinsZ-1;
  if(bin1 == bin2) {
#ifdef G4ANALYSIS_USE
    if(tree) {
      histo[3]->fill(z0,edep);
      if(r1 < stepR) {
        G4double w = edep;
        if(r2 > stepR) w *= (stepR - r1)/(r2 - r1);
        histo[4]->fill(z0,w);
      }
    }
#endif
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
    G4double zz2 = (bin1+1)*stepZ;
    G4double rr1 = r1;
    G4double dz  = z2 - z1;
    G4double dr  = r2 - r1;
    G4double rr2 = r1 + dr*(zz2-zz1)/dz;
    for(bin=bin1; bin<=bin2; bin++) {
#ifdef G4ANALYSIS_USE
      if(tree) {
	G4double de = edep*(zz2 - zz1)/dz;
	G4double zf = (zz1+zz2)*0.5;
	histo[3]->fill(zf,de);
	if(rr1 < stepR) {
	  G4double w = de;
	  if(rr2 > stepR) w *= (stepR - rr1)/(rr2 - rr1);
	  histo[4]->fill(zf,w);
	}
      }
#endif
      zz1 = zz2;
      zz2 = std::min(z2, zz1+stepZ);
      rr1 = rr2;
      rr2 = rr1 + dr*(zz2 - zz1)/dz;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

