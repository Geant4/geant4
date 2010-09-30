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

#include "HistoITEPTest47.hh"
#include "G4VParticleChange.hh"
#include "G4UnitsTable.hh"

HistoITEPTest47::HistoITEPTest47(std::string namePart, std::string nameMat, 
				 G4double momentum, std::string nameGen) :
  HistoTest47(namePart, nameMat, momentum, nameGen) {

  energies.push_back(0.01);
  energies.push_back(0.03);
  energies.push_back(0.05);
  energies.push_back(0.07);
  energies.push_back(0.09);
  energies.push_back(0.11);
  energies.push_back(0.13);
  energies.push_back(0.15);
  energies.push_back(0.17);
  energies.push_back(0.19);
  energies.push_back(0.21);
  energies.push_back(0.23);
  energies.push_back(0.25);
  angles.push_back(10.1*deg);
  angles.push_back(15.0*deg);
  angles.push_back(19.8*deg);
  angles.push_back(24.8*deg);
  angles.push_back(29.5*deg);
  angles.push_back(34.6*deg);
  angles.push_back(39.6*deg);
  angles.push_back(44.3*deg);
  angles.push_back(49.3*deg);
  angles.push_back(54.2*deg);
  angles.push_back(59.1*deg);
  angles.push_back(64.1*deg);
  angles.push_back(69.1*deg);
  angles.push_back(74.1*deg);
  angles.push_back(79.1*deg);
  angles.push_back(84.1*deg);
  angles.push_back(89.0*deg);
  angles.push_back(98.9*deg);
  angles.push_back(108.9*deg);
  angles.push_back(119.0*deg);
  angles.push_back(129.1*deg);
  angles.push_back(139.1*deg);
  angles.push_back(149.3*deg);
  angles.push_back(159.6*deg);
  angles.push_back(161.4*deg);
  angles.push_back(165.5*deg);
  angles.push_back(169.5*deg);
  angles.push_back(173.5*deg);
  angles.push_back(177.0*deg);

  dtheta = 4.0*deg;
  de     = 0.02;

  for (unsigned int ii=0; ii<angles.size(); ii++) {

    double cth1 = std::cos(std::min(angles[ii]+dtheta,180*deg));
    double cth2 = std::cos(std::max(angles[ii]-dtheta,0.0*deg));
    dcth.push_back(std::abs(cth1-cth2));
    cthmin.push_back(std::min(cth1,cth2));
    cthmax.push_back(std::max(cth1,cth2));
  }

  for (unsigned int ii=0; ii<energies.size(); ii++) {
    double en = energies[ii];
    emin.push_back(en-0.5*de);
    emax.push_back(en+0.5*de);
  }

  G4cout << "HistoITEPTest47:: Initialized with " << cthmin.size() << " theta bins and " << emin.size() << " energy bins" << G4endl;
  for (unsigned int ii=0; ii<cthmin.size(); ii++) 
    G4cout << "HistoITEPTest47::Bin " << ii << " theta " << angles[ii]/deg << " cos(theta) = " << cthmin[ii] << ":" << cthmax[ii] << " dcostheta " << dcth[ii] << G4endl;
  for (unsigned int ii=0; ii<emin.size(); ii++) 
    G4cout << "HistoITEPTest47::Bin " << ii << " Energy = " << emin[ii] << ":" << emax[ii] << G4endl;

}

HistoITEPTest47::~HistoITEPTest47() {}

void HistoITEPTest47::fill(G4VParticleChange* aChange, G4LorentzVector pinit) {

  if (unInitialized) initialize();

  G4LorentzVector labv(pinit), fm;
  G4int n = aChange->GetNumberOfSecondaries();
  const G4DynamicParticle* sec = 0;
  G4ParticleDefinition* pd;
  G4ThreeVector  mom;
  for (G4int j=0; j<n; j++) {
    sec = aChange->GetSecondary(j)->GetDynamicParticle();
    pd  = sec->GetDefinition();
    mom = sec->GetMomentumDirection();
    G4double ke   = (sec->GetKineticEnergy())/GeV;
    if (ke < 0.0) ke = 0.0;
    G4double m     = (pd->GetPDGMass())/GeV;
    G4double p     = std::sqrt(ke*(ke + 2.0*m));
    mom           *= p;
    fm             = G4LorentzVector(mom, ke+m);
    labv          -= fm;
    G4int    type  = particleType(pd);
    if (type >= 0 && type <= 5) {
      G4double theta = mom.theta();
      G4double cth   = std::cos(theta);
      G4double wt    = 1.0/p;
      for (unsigned int ii=0; ii<angles.size(); ii++) {
	if (cth > cthmin[ii] && cth <= cthmax[ii]) {
          switch (type) {
          case 0:
	    hiKE11[ii]->Fill(ke);
	    hiKE12[ii]->Fill(ke,wt);
            break;
          case 1:
	    hiKE21[ii]->Fill(ke);
	    hiKE22[ii]->Fill(ke,wt);
            break;
          case 2:
	    hiKE31[ii]->Fill(ke);
	    hiKE32[ii]->Fill(ke,wt);
            break;
          case 3:
	    hiKE41[ii]->Fill(ke);
	    hiKE42[ii]->Fill(ke,wt);
            break;
          case 4:
	    hiKE51[ii]->Fill(ke);
	    hiKE52[ii]->Fill(ke,wt);
            break;
          case 5:
	    hiKE61[ii]->Fill(ke);
	    hiKE62[ii]->Fill(ke,wt);
	    break;
	  }
	}
      }
      if (type <= 1) {
	for (unsigned int ii=0; ii<energies.size(); ii++) {
	  if (ke > emin[ii] && ke <= emax[ii]) {
	    if (type == 0) {
	      hiCT11[ii]->Fill(cth);
	      hiCT12[ii]->Fill(cth,wt);
	    } else {
	      hiCT21[ii]->Fill(cth);
	      hiCT22[ii]->Fill(cth,wt);
	    }
	  }
	}
      }
    }
  }

  epTest.fill(aChange,pinit);
}

void HistoITEPTest47::write(G4double cross_sec, G4int nevt) {

  char name[100], title[100];
  G4double xbin, scale;
  std::vector<TH1F*> hiKE10, hiKE20, hiKE30, hiKE40, hiKE50, hiKE60;
  for (unsigned int ii=0; ii<angles.size(); ii++) {
    xbin = hiKE11[ii]->GetBinWidth(1);
    sprintf (title, "Kinetic Energy of p (GeV)");
    hiKE11[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiKE11[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiKE12[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dcth[ii]);
    sprintf (title, "Kinetic Energy of p (GeV)");
    hiKE12[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiKE12[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "KEproton0%s%5.1f", tag1Name,  angles[ii]/deg);
    hiKE10.push_back((TH1F*)hiKE12[ii]->Clone());
    hiKE10[ii]->SetName(name); 
    hiKE10[ii]->Scale(scale);
    hiKE10[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiKE21[ii]->GetBinWidth(1);
    sprintf (title, "Kinetic Energy of n (GeV)");
    hiKE21[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiKE21[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiKE22[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dcth[ii]);
    sprintf (title, "Kinetic Energy of n (GeV)");
    hiKE22[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiKE22[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "KEneutron0%s%5.1f", tag1Name,  angles[ii]/deg);
    hiKE20.push_back((TH1F*)hiKE22[ii]->Clone());
    hiKE20[ii]->SetName(name); 
    hiKE20[ii]->Scale(scale);
    hiKE20[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiKE31[ii]->GetBinWidth(1);
    sprintf (title, "Kinetic Energy of #pi+ (GeV)");
    hiKE31[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiKE31[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiKE32[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dcth[ii]);
    sprintf (title, "Kinetic Energy of #pi+ (GeV)");
    hiKE32[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiKE32[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "KEpiplus0%s%5.1f", tag1Name,  angles[ii]/deg);
    hiKE30.push_back((TH1F*)hiKE32[ii]->Clone());
    hiKE30[ii]->SetName(name); 
    hiKE30[ii]->Scale(scale);
    hiKE30[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiKE41[ii]->GetBinWidth(1);
    sprintf (title, "Kinetic Energy of #pi- (GeV)");
    hiKE41[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiKE41[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiKE42[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dcth[ii]);
    sprintf (title, "Kinetic Energy of #pi- (GeV)");
    hiKE42[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiKE42[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "KEpiminus0%s%5.1f", tag1Name,  angles[ii]/deg);
    hiKE40.push_back((TH1F*)hiKE42[ii]->Clone());
    hiKE40[ii]->SetName(name); 
    hiKE40[ii]->Scale(scale);
    hiKE40[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiKE51[ii]->GetBinWidth(1);
    sprintf (title, "Kinetic Energy of K+ (GeV)");
    hiKE51[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiKE51[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiKE52[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dcth[ii]);
    sprintf (title, "Kinetic Energy of K+ (GeV)");
    hiKE52[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiKE52[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "KEkplus0%s%5.1f", tag1Name,  angles[ii]/deg);
    hiKE50.push_back((TH1F*)hiKE52[ii]->Clone());
    hiKE50[ii]->SetName(name); 
    hiKE50[ii]->Scale(scale);
    hiKE50[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiKE61[ii]->GetBinWidth(1);
    sprintf (title, "Kinetic Energy of K- (GeV)");
    hiKE61[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiKE61[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiKE62[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dcth[ii]);
    sprintf (title, "Kinetic Energy of K- (GeV)");
    hiKE62[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiKE62[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "KEkminus0%s%5.1f", tag1Name,  angles[ii]/deg);
    hiKE60.push_back((TH1F*)hiKE62[ii]->Clone());
    hiKE60[ii]->SetName(name); 
    hiKE60[ii]->Scale(scale);
    hiKE60[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");
  }

  std::vector<TH1F*> hiCT10, hiCT20;
  for (unsigned int ii=0; ii<energies.size(); ii++) {
    xbin = hiCT11[ii]->GetBinWidth(1);
    sprintf (title, "Events/%6.3f", xbin);
    hiCT11[ii]->GetXaxis()->SetTitle("cos (#theta)");
    hiCT11[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiCT12[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*de);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f", xbin);
    hiCT12[ii]->GetXaxis()->SetTitle("cos (#theta)");
    hiCT12[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "CTproton0%s%4.2f", tag1Name, energies[ii]);
    hiCT10.push_back((TH1F*)hiCT12[ii]->Clone());
    hiCT10[ii]->SetName(name); 
    hiCT10[ii]->Scale(scale);
    hiCT10[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiCT12[ii]->GetBinWidth(1);
    sprintf (title, "Events/%6.3f", xbin);
    hiCT12[ii]->GetXaxis()->SetTitle("cos (#theta)");
    hiCT12[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiCT22[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*de);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f", xbin);
    hiCT22[ii]->GetXaxis()->SetTitle("cos (#theta)");
    hiCT22[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "CTneutron0%s%4.2f", tag1Name, energies[ii]);
    hiCT20.push_back((TH1F*)hiCT22[ii]->Clone());
    hiCT20[ii]->SetName(name); 
    hiCT20[ii]->Scale(scale);
    hiCT20[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");
  }

  TFile f(fileName.c_str(), "recreate");
  for (unsigned int ii=0; ii<angles.size(); ii++) {
    hiKE11[ii]->Write(); hiKE10[ii]->Write(); hiKE12[ii]->Write();
    hiKE21[ii]->Write(); hiKE20[ii]->Write(); hiKE22[ii]->Write();
    hiKE31[ii]->Write(); hiKE30[ii]->Write(); hiKE32[ii]->Write();
    hiKE41[ii]->Write(); hiKE40[ii]->Write(); hiKE42[ii]->Write();
    hiKE51[ii]->Write(); hiKE50[ii]->Write(); hiKE52[ii]->Write();
    hiKE61[ii]->Write(); hiKE60[ii]->Write(); hiKE62[ii]->Write();
  }
  for (unsigned int ii=0; ii<energies.size(); ii++) {
    hiCT11[ii]->Write(); hiCT10[ii]->Write(); hiCT12[ii]->Write();
    hiCT21[ii]->Write(); hiCT20[ii]->Write(); hiCT22[ii]->Write();
  }
  epTest.write();
  f.Close();
}

void HistoITEPTest47::initialize() {

  unInitialized = false;
  G4cout << "HistoITEPTest47::initialize invoked" << G4endl;
  fileName = particle + target + generator;
  std::ostringstream tmp;
  char nams[4];
  sprintf (nams, "%4.2f", energy);
  tmp << nams << "GeV";
  if ( jobID > -1 )     tmp << "-" << jobID;
  if ( clusterID > -1 ) tmp << "-" << clusterID;
  tmp << ".root";
  fileName += tmp.str();
  sprintf (tag1Name, "%s%s%s%4.2fGeV", particle.c_str(), target.c_str(),
	   generator.c_str(), energy); 
  sprintf (tag2Name, "%s+%s", particle.c_str(), target.c_str());
  sprintf (tag3Name, "at %4.2f GeV (%s)", energy, generator.c_str());
  G4cout << "HistoITEPTest47::fileName:" << fileName.c_str() 
	 << " Tag1:" << tag1Name << " Tag2: " << tag2Name << " Tag3: " 
	 << tag3Name << G4endl;

  book();
  epTest.initialize(particle,target,energy,generator);
}

void HistoITEPTest47::book() {

  char name[100], title[160];
  hiKE11.clear(); hiKE12.clear(); hiKE21.clear(); hiKE22.clear();
  hiKE31.clear(); hiKE32.clear(); hiKE41.clear(); hiKE42.clear();
  hiKE51.clear(); hiKE52.clear(); hiKE61.clear(); hiKE62.clear();
  for (unsigned int ii=0; ii<angles.size(); ii++) {
    double ang = angles[ii]/deg;
    sprintf (name, "KEproton1%s%5.1f", tag1Name,  ang);
    sprintf (title, "%s to p %s (#theta = %8.2f)", tag2Name, tag3Name, ang);
    hiKE11.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEproton2%s%5.1f", tag1Name,  ang);
    hiKE12.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEneutron1%s%5.1f", tag1Name,  ang);
    sprintf (title, "%s to n %s (#theta = %8.2f)", tag2Name, tag3Name, ang);
    hiKE21.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEneutron2%s%5.1f", tag1Name,  ang);
    hiKE22.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEpiplus1%s%5.1f", tag1Name,  ang);
    sprintf (title, "%s to #pi+ %s (#theta = %8.2f)", tag2Name, tag3Name, ang);
    hiKE31.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEpiplu2%s%5.1f", tag1Name,  ang);
    hiKE32.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEpiminus1%s%5.1f", tag1Name,  ang);
    sprintf (title, "%s to #pi- %s (#theta = %8.2f)", tag2Name, tag3Name, ang);
    hiKE41.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEpiminus2%s%5.1f", tag1Name,  ang);
    hiKE42.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEkplus1%s%5.1f", tag1Name,  ang);
    sprintf (title, "%s to K+ %s (#theta = %8.2f)", tag2Name, tag3Name, ang);
    hiKE51.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEkplu2%s%5.1f", tag1Name,  ang);
    hiKE52.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEkminus1%s%5.1f", tag1Name,  ang);
    sprintf (title, "%s to K- %s (#theta = %8.2f)", tag2Name, tag3Name, ang);
    hiKE61.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KElminus2%s%5.1f", tag1Name,  ang);
    hiKE62.push_back(new TH1F (name, title, 800, 0., 8.));
  }

  hiCT11.clear (); hiCT12.clear(); hiCT21.clear(); hiCT22.clear();
  for (unsigned int ii=0; ii<energies.size(); ii++) {
    double en = energies[ii];
    emin.push_back(en-0.5*de);
    emax.push_back(en+0.5*de);
    sprintf (name, "CTproton1%s%4.2f", tag1Name, energies[ii]);
    sprintf (title, "%s to p at %s (KE = %6.2f GeV)", tag2Name, tag3Name,
	     energies[ii]);
    hiCT11.push_back(new TH1F (name, title, 80, -1., 1.));
    sprintf (name, "CTproton2%s%4.2f", tag1Name, energies[ii]);
    hiCT12.push_back(new TH1F (name, title, 80, -1., 1.));
    sprintf (name, "CTneutron1%s%4.2f", tag1Name, energies[ii]);
    sprintf (title, "%s to n at %s (KE = %6.2f GeV)", tag2Name, tag3Name,
	     energies[ii]);
    hiCT21.push_back(new TH1F (name, title, 80, -1., 1.));
    sprintf (name, "CTneutron2%s%4.2f", tag1Name, energies[ii]);
    hiCT22.push_back(new TH1F (name, title, 80, -1., 1.));
  }
}
