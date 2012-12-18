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
#include "globals.hh"
#include "G4ios.hh"

#include "HistoMIPSTest47.hh"
#include "G4VParticleChange.hh"
#include "G4UnitsTable.hh"
#include "CLHEP/Units/PhysicalConstants.h"

HistoMIPSTest47::HistoMIPSTest47(std::string namePart, std::string nameMat, 
				 G4double momentum, std::string nameGen) :
  HistoTest47(namePart, nameMat, momentum, nameGen) {


  xminCalorimeter = -49.5*CLHEP::cm;
  xmaxCalorimeter =  49.5*CLHEP::cm;
  yminCalorimeter = -49.0*CLHEP::cm;
  ymaxCalorimeter =  49.0*CLHEP::cm;
  rmaxCalorimeter =  45.0*CLHEP::cm;
  zCalorimeter    = 2596.0*CLHEP::cm;
  sterad          = 0.0013;

  G4cout << "HistoMIPSTest47:: Initialized with Calorimeter position "
	 << zCalorimeter << " Ranges (rmax:xmin:xmax:ymin:ymax) " 
	 << rmaxCalorimeter << ":" << xminCalorimeter << ":" << xmaxCalorimeter
	 << ":" << yminCalorimeter << ":" << ymaxCalorimeter << " sterad " 
	 << sterad << G4endl;

}

HistoMIPSTest47::~HistoMIPSTest47() {}

void HistoMIPSTest47::fill(G4VParticleChange* aChange, G4LorentzVector& pinit,
			   G4ThreeVector& aPosition, G4LorentzVector& labp) {

  if (unInitialized) initialize();

  G4LorentzVector labv(pinit), fm;
  const G4DynamicParticle* sec = 0;
  G4ParticleDefinition* pd;
  G4ThreeVector  mom;
  G4int          n       = aChange->GetNumberOfSecondaries();
  G4double       distz   = zCalorimeter - aPosition.z();
  G4ThreeVector  boostCM = labp.findBoostToCM();
  G4double       roots   = labp.mag();
  if (debug) G4cout << "HistoMIPSTest47::fill CM system " << labv << " | "
		    << labp << " Boost " << boostCM << " roots " << roots
		    << " N " << n  << G4endl;
  for (G4int j=0; j<n; j++) {
    sec = aChange->GetSecondary(j)->GetDynamicParticle();
    pd  = sec->GetDefinition();
    mom = sec->GetMomentumDirection();
    G4double ke   = (sec->GetKineticEnergy())/CLHEP::GeV;
    if (ke < 0.0) ke = 0.0;
    G4double m     = (pd->GetPDGMass())/CLHEP::GeV;
    G4double p     = std::sqrt(ke*(ke + 2.0*m));
    G4double pmax  = std::sqrt(0.25*roots*roots-m*m);
    mom           *= p;
    fm             = G4LorentzVector(mom, ke+m);
    labv          -= fm;
    G4int    type  = particleType(pd);
    G4bool   acc1=false, acc2=false;
    if (debug) G4cout << "Particle " << j << " type " << type << " p " << mom
		      << " p|m|pmax " << p << "|" << m << "|" << pmax <<G4endl;
    if (mom.z() > 0 && p > eMin) {
      G4double xCal = (distz*mom.x()/mom.z()) + aPosition.x();
      G4double yCal = (distz*mom.y()/mom.z()) + aPosition.y();
      G4double rCal = std::sqrt(xCal*xCal+yCal*yCal);
      if (rCal < rmaxCalorimeter) acc1 = true;
      if (xCal > xminCalorimeter && xCal < xmaxCalorimeter &&
	  yCal > yminCalorimeter && yCal < ymaxCalorimeter) acc2 = true;
      if (debug) G4cout << "x/y/r at Calorimeter " << xCal << ":" << yCal
			<< ":" << rCal << " Accept " << acc1 << ":" << acc2
			<< G4endl;
    }
    if (type >= 0 && type <= 7) {
      G4LorentzVector fmc = fm.boost(boostCM.x(),boostCM.y(),boostCM.z());
      G4double        xf  = fmc.z()/pmax;
      G4double        wt  = (fmc.e()*fmc.z())/(fmc.rho()*fmc.rho()*fmc.rho()*pmax);
      if (debug) G4cout << "Type " << type << " CM Momentum " << fmc << " p|pt " 
			<< fmc.rho() << "|" << fmc.perp() << "|" << pmax 
			<< " xf " << xf << " wt " << wt << " Acc: "
			<< acc1 << ":" << acc2 << G4endl;
      hiPL1[type]->Fill(p);
      hiXF1[type]->Fill(xf);
      hiXW1[type]->Fill(xf,wt);
      if (acc1) {
	hiPL2[type]->Fill(p);
	hiXF2[type]->Fill(xf);
	hiXW2[type]->Fill(xf,wt);
      }
      if (acc2) {
	hiPL3[type]->Fill(p);
	hiXF3[type]->Fill(xf);
	hiXW3[type]->Fill(xf,wt);
      }
    }
  }

  epTest.fill(aChange,pinit,aPosition);
}

void HistoMIPSTest47::write(G4double cross_sec, G4int nevt) {

  char name[100], title[100];
  G4double xbin, scale;
  std::vector<TH1F*> hiPL0, hiXF0, hiXW0;
  std::vector<TH1F*> hiPLx, hiXFx, hiXWx;
  std::string psymbs[8]={"p","n","#pi^{+}","#pi^{-}","K^{+}","K^{-}",
			 "pbar","#pi^{0}"};
  std::string particles[8]={"proton", "neutron", "piplus", "piminus",
			    "Kplus", "Kminus", "pbar", "pizero"};
  for (unsigned int ii=0; ii<8; ii++) {
    xbin = hiPL1[ii]->GetBinWidth(1);
    sprintf (title, "Laboratory momentum of %s (GeV/c)", psymbs[ii].c_str());
    hiPL1[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.2f GeV/c", xbin);
    hiPL1[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiPL2[ii]->GetBinWidth(1);
    sprintf (title, "Laboratory momentum of %s (GeV/c)", psymbs[ii].c_str());
    hiPL2[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.2f GeV/c", xbin);
    hiPL2[ii]->GetYaxis()->SetTitle(title);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin);
    sprintf (name, "PL%s0%s", particles[ii].c_str(), tag1Name);
    hiPL0.push_back((TH1F*)hiPL2[ii]->Clone());
    hiPL0[ii]->SetName(name); 
    hiPL0[ii]->Scale(scale);
    hiPL0[ii]->GetYaxis()->SetTitle("#frac{d#sigma}{dp} (mb/GeV)");
    xbin = hiPL3[ii]->GetBinWidth(1);
    sprintf (title, "Laboratory momentum of %s (GeV/c)", psymbs[ii].c_str());
    hiPL3[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.2f GeV/c", xbin);
    hiPL3[ii]->GetYaxis()->SetTitle(title);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin);
    sprintf (name, "PL%sx%s", particles[ii].c_str(), tag1Name);
    hiPLx.push_back((TH1F*)hiPL3[ii]->Clone());
    hiPLx[ii]->SetName(name); 
    hiPLx[ii]->Scale(scale);
    hiPLx[ii]->GetYaxis()->SetTitle("#frac{d#sigma}{dp} (mb/GeV)");

    xbin = hiXF1[ii]->GetBinWidth(1);
    sprintf (title, "x_{F} of %s", psymbs[ii].c_str());
    hiXF1[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%5.2f", xbin);
    hiXF1[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiXF2[ii]->GetBinWidth(1);
    sprintf (title, "x_{F} of %s", psymbs[ii].c_str());
    hiXF2[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%5.2f", xbin);
    hiXF2[ii]->GetYaxis()->SetTitle(title);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin);
    sprintf (name, "XF%s0%s", particles[ii].c_str(), tag1Name);
    hiXF0.push_back((TH1F*)hiXF2[ii]->Clone());
    hiXF0[ii]->SetName(name); 
    hiXF0[ii]->Scale(scale);
    hiXF0[ii]->GetYaxis()->SetTitle(title);
    xbin = hiXF3[ii]->GetBinWidth(1);
    sprintf (title, "x_{F} of %s", psymbs[ii].c_str());
    hiXF3[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%5.2f", xbin);
    hiXF3[ii]->GetYaxis()->SetTitle(title);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin);
    sprintf (name, "XF%sx%s", particles[ii].c_str(), tag1Name);
    hiXFx.push_back((TH1F*)hiXF3[ii]->Clone());
    hiXFx[ii]->SetName(name); 
    hiXFx[ii]->Scale(scale);
    hiXFx[ii]->GetYaxis()->SetTitle(title);

    xbin = hiXW1[ii]->GetBinWidth(1);
    sprintf (title, "x_{F} of %s", psymbs[ii].c_str());
    hiXW1[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{E}{p^{2}})/%5.2f", xbin);
    hiXW1[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiXW2[ii]->GetBinWidth(1);
    sprintf (title, "x_{F} of %s", psymbs[ii].c_str());
    hiXW2[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{E}{p^{2}})/%5.2f", xbin);
    hiXW2[ii]->GetYaxis()->SetTitle(title);
    scale = (4.*CLHEP::pi*cross_sec)/(((double)(std::max(nevt,1)))*xbin);
    sprintf (name, "XW%s0%s", particles[ii].c_str(), tag1Name);
    hiXW0.push_back((TH1F*)hiXW2[ii]->Clone());
    hiXW0[ii]->SetName(name); 
    hiXW0[ii]->Scale(scale);
    hiXW0[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");
    xbin  = hiXW3[ii]->GetBinWidth(1);
    sprintf (title, "x_{F} of %s", psymbs[ii].c_str());
    hiXW3[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{E}{p^{2}})/%5.2f", xbin);
    hiXW3[ii]->GetYaxis()->SetTitle(title);
    scale = (4.*CLHEP::pi*cross_sec)/(((double)(std::max(nevt,1)))*xbin);
    sprintf (name, "XW%sx%s", particles[ii].c_str(), tag1Name);
    hiXWx.push_back((TH1F*)hiXW3[ii]->Clone());
    hiXWx[ii]->SetName(name); 
    hiXWx[ii]->Scale(scale);
    hiXWx[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");
  }

  TFile f(fileName.c_str(), "recreate");
  for (unsigned int ii=0; ii<8; ii++) {
    hiPL0[ii]->Write(); hiPL1[ii]->Write(); hiPL2[ii]->Write(); hiPL3[ii]->Write(); hiPLx[ii]->Write(); 
    hiXF0[ii]->Write(); hiXF1[ii]->Write(); hiXF2[ii]->Write(); hiXF3[ii]->Write(); hiXFx[ii]->Write();
    hiXW0[ii]->Write(); hiXW1[ii]->Write(); hiXW2[ii]->Write(); hiXW3[ii]->Write(); hiXWx[ii]->Write(); 
  }
  epTest.write();
  f.Close();
}

void HistoMIPSTest47::initialize() {

  unInitialized = false;
  if (energy < 60.0) {
    eMin     = 12.0;
  } else if (energy > 100.0) {
    eMin     = 20.0;
  } else {
    eMin     = 18.0;
  }
  G4cout << "HistoMIPSTest47::initialize invoked: E(eMin) "  << energy << "(" 
	 << eMin << ")" << G4endl;
  fileName = particle + target + generator;
  std::ostringstream tmp;
  char nams[4];
  int iene = (int)(energy);
  if (iene >= 100) sprintf (nams, "%3d.",  iene);
  else             sprintf (nams, "%2d.0", iene);
  tmp << nams << "GeV";
  if ( jobID > -1 )     tmp << "-" << jobID;
  if ( clusterID > -1 ) tmp << "-" << clusterID;
  tmp << ".root";
  fileName += tmp.str();
  sprintf (tag1Name, "%s%s%s%sGeV", particle.c_str(), target.c_str(), 
	   generator.c_str(), nams); 
  sprintf (tag2Name, "%s+%s", particle.c_str(), target.c_str());
  sprintf (tag3Name, "at %s GeV (%s)", nams, generator.c_str());
  G4cout << "HistoMIPSTest47::fileName:" << fileName.c_str() 
	 << " Tag1:" << tag1Name << " Tag2: " << tag2Name << " Tag3: " 
	 << tag3Name << " eMin: " << eMin << G4endl;

  book();
  epTest.initialize(particle,target,energy,generator);
  epTest.setDebug(debug);
}

void HistoMIPSTest47::book() {

  char name[100], title[160];
  hiPL1.clear(); hiPL2.clear(); hiPL3.clear();
  hiXF1.clear(); hiXF2.clear(); hiXF3.clear();
  hiXW1.clear(); hiXW2.clear(); hiXW3.clear();
  std::string particles[8]={"proton", "neutron", "piplus", "piminus",
			    "Kplus", "Kminus", "pbar", "pizero"};
  for (unsigned int ii=0; ii<8; ii++) {
    sprintf (name, "PL%s1%s", particles[ii].c_str(), tag1Name);
    sprintf (title, "%s to %s %s", tag2Name, particles[ii].c_str(), tag3Name);
    hiPL1.push_back(new TH1F (name, title, 150, 0., 150.));
    hiPL1[ii]->Sumw2();
    sprintf (name, "PL%s2%s", particles[ii].c_str(), tag1Name);
    hiPL2.push_back(new TH1F (name, title, 150, 0., 150.));
    hiPL2[ii]->Sumw2();
    sprintf (name, "PL%s3%s", particles[ii].c_str(), tag1Name);
    hiPL3.push_back(new TH1F (name, title, 150, 0., 150.));
    hiPL3[ii]->Sumw2();
    sprintf (name, "XF%s1%s", particles[ii].c_str(), tag1Name);
    hiXF1.push_back(new TH1F (name, title, 40, -1.0, 1.0));
    hiXF1[ii]->Sumw2();
    sprintf (name, "XF%s2%s", particles[ii].c_str(), tag1Name);
    hiXF2.push_back(new TH1F (name, title, 40, -1.0, 1.0));
    hiXF2[ii]->Sumw2();
    sprintf (name, "XF%s3%s", particles[ii].c_str(), tag1Name);
    hiXF3.push_back(new TH1F (name, title, 40, -1.0, 1.0));
    hiXF3[ii]->Sumw2();
    sprintf (name, "XW%s1%s", particles[ii].c_str(), tag1Name);
    hiXW1.push_back(new TH1F (name, title, 40, -1.0, 1.0));
    hiXW1[ii]->Sumw2();
    sprintf (name, "XW%s2%s", particles[ii].c_str(), tag1Name);
    hiXW2.push_back(new TH1F (name, title, 40, -1.0, 1.0));
    hiXW2[ii]->Sumw2();
    sprintf (name, "XW%s3%s", particles[ii].c_str(), tag1Name);
    hiXW3.push_back(new TH1F (name, title, 40, -1.0, 1.0));
    hiXW3[ii]->Sumw2();
  }
}
