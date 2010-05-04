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

    double cth1 = cos(std::min(angles[ii]+dtheta,180*deg));
    double cth2 = cos(std::max(angles[ii]-dtheta,0.0*deg));
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
    G4double p     = sqrt(ke*(ke + 2.0*m));
    mom           *= p;
    fm             = G4LorentzVector(mom, ke+m);
    labv          -= fm;
    G4int    type  = particleType(pd);
    if (type >= 0 && type <= 1) {
      G4double theta = mom.theta();
      G4double cth   = cos(theta);
      G4double wt    = 1.0/p;
      for (unsigned int ii=0; ii<angles.size(); ii++) {
	if (cth > cthmin[ii] && cth <= cthmax[ii]) {
	  if (type == 0) {
	    hiKE11[ii]->Fill(ke);
	    hiKE12[ii]->Fill(ke,wt);
	  } else {
	    hiKE21[ii]->Fill(ke);
	    hiKE22[ii]->Fill(ke,wt);
	  }
	}
      }
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

  epTest.fill(aChange,pinit);
}

void HistoITEPTest47::write(G4double cross_sec, G4int nevt) {

  char name[100], title[100];
  G4double xbin, scale;
  std::vector<TH1F*> hiKE10, hiKE20;
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
  for (unsigned int ii=0; ii<angles.size(); ii++) {
    sprintf (name, "KEproton1%s%5.1f", tag1Name,  angles[ii]/deg);
    sprintf (title, "%s to p %s (#theta = %8.2f)", tag2Name, tag3Name,
	     angles[ii]/deg);
    hiKE11.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEproton2%s%5.1f", tag1Name,  angles[ii]/deg);
    hiKE12.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEneutron1%s%5.1f", tag1Name,  angles[ii]/deg);
    sprintf (title, "%s to n %s (#theta = %8.2f)", tag2Name, tag3Name,
	     angles[ii]/deg);
    hiKE21.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "KEneutron2%s%5.1f", tag1Name,  angles[ii]/deg);
    hiKE22.push_back(new TH1F (name, title, 800, 0., 8.));
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
