#include "globals.hh"
#include "G4ios.hh"

#include "HistoBNLTest47.hh"
#include "G4VParticleChange.hh"
#include "G4UnitsTable.hh"

HistoBNLTest47::HistoBNLTest47(std::string namePart, std::string nameMat, 
			       G4double momentum, std::string nameGen) :
  HistoTest47(namePart, nameMat, momentum, nameGen) {

  rapidities.push_back(0.60);
  rapidities.push_back(0.80);
  rapidities.push_back(1.00);
  rapidities.push_back(1.20);
  rapidities.push_back(1.40);
  rapidities.push_back(1.60);
  rapidities.push_back(1.80);
  rapidities.push_back(2.00);
  rapidities.push_back(2.20);
  rapidities.push_back(2.40);
  rapidities.push_back(2.60);
  rapidities.push_back(2.80);
  rapidities.push_back(3.00);

  for (unsigned int ii=0; ii<rapidities.size()-1; ii++) {
    ymin.push_back(rapidities[ii]);
    ymax.push_back(rapidities[ii+1]);
  }

  G4cout << "HistoBNLTest47:: Initialized with " << ymin.size() << " rapidity bins" << G4endl;
  for (unsigned int ii=0; ii<ymin.size(); ii++) 
    G4cout << "HistoBNLTest47::Bin " << ii << " rapidity = " << ymin[ii] << ":" << ymax[ii] << G4endl;

}

HistoBNLTest47::~HistoBNLTest47() {}

void HistoBNLTest47::fill(G4VParticleChange* aChange, G4LorentzVector pinit) {

  if (unInitialized) initialize();

  G4LorentzVector labv(pinit), fm;
  G4int n = aChange->GetNumberOfSecondaries();
  const G4DynamicParticle* sec = 0;
  G4ThreeVector  mom;
  G4ParticleDefinition* pd;
  for (G4int j=0; j<n; j++) {
    sec = aChange->GetSecondary(j)->GetDynamicParticle();
    pd  = sec->GetDefinition();
    mom = sec->GetMomentumDirection();
    G4int    type  = particleType(pd);
    G4double ke    = (sec->GetKineticEnergy())/GeV;
    if (ke < 0.0) ke = 0.0;
    G4double m     = (pd->GetPDGMass())/GeV;
    G4double p     = sqrt(ke*(ke + 2.0*m));
    G4double ee    = ke + m;
    mom           *= p;
    fm             = G4LorentzVector(mom, ee);
    labv          -= fm;
    if (type >= 0 && type <= 5) {
      G4double theta = mom.theta();
      G4double pt    = p*sin(theta);
      G4double pl    = p*cos(theta);
      G4double mt    = sqrt (pt*pt + m*m);
      G4double mtp   = (mt - std::abs(m));
      G4double yv    = 0.5*log((ee+pl)/(ee-pl));
      G4double wt    = 1./mt;
      for (unsigned int ii=0; ii<ymin.size(); ii++) {
	if (yv > ymin[ii] && yv <= ymax[ii]) {
	  switch (type) {
	  case 0:
	    hiMT11[ii]->Fill(mtp);
	    hiMT12[ii]->Fill(mtp,wt);
	    break;
	  case 1:
	    hiMT21[ii]->Fill(mtp);
	    hiMT22[ii]->Fill(mtp,wt);
	    break;
	  case 2:
	    hiMT31[ii]->Fill(mtp);
	    hiMT32[ii]->Fill(mtp,wt);
	    break;
	  case 3:
	    hiMT41[ii]->Fill(mtp);
	    hiMT42[ii]->Fill(mtp,wt);
	    break;
	  case 4:
	    hiMT51[ii]->Fill(mtp);
	    hiMT52[ii]->Fill(mtp,wt);
	    break;
	  case 5:
	    hiMT61[ii]->Fill(mtp);
	    hiMT62[ii]->Fill(mtp,wt);
	    break;
	  }
	}
      }
    }
  }

  epTest.fill(aChange,pinit);
}

void HistoBNLTest47::write(G4double cross_sec, G4int nevt) {

  char name[100], title[100];
  G4double xbin, dy, yv, scale;
  std::vector<TH1F*> hiMT10, hiMT20, hiMT30, hiMT40, hiMT50, hiMT60;
  for (unsigned int ii=0; ii<ymin.size(); ii++) {
    dy    = (ymax[ii]-ymin[ii]);
    yv    = 0.5*(ymin[ii]+ymax[ii]);

    xbin  = hiMT11[ii]->GetBinWidth(1);
    sprintf (title, "Reduced Transverse Mass of p (GeV)");
    hiMT11[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiMT11[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiMT12[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dy);
    sprintf (title, "Reduced Transverse Mass of p (GeV)");
    hiMT12[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiMT12[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "MTproton0%s%4.2f", tag1Name,  yv);
    hiMT10.push_back((TH1F*)hiMT12[ii]->Clone());
    hiMT10[ii]->SetName(name); 
    hiMT10[ii]->Scale(scale);
    hiMT10[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin  = hiMT21[ii]->GetBinWidth(1);
    sprintf (title, "Reduced Transverse Mass of n (GeV)");
    hiMT21[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiMT21[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiMT22[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dy);
    sprintf (title, "Reduced Transverse Mass of n (GeV)");
    hiMT22[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiMT22[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "MTneutron0%s%4.2f", tag1Name,  yv);
    hiMT20.push_back((TH1F*)hiMT22[ii]->Clone());
    hiMT20[ii]->SetName(name); 
    hiMT20[ii]->Scale(scale);
    hiMT20[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin  = hiMT31[ii]->GetBinWidth(1);
    sprintf (title, "Reduced Transverse Mass of #pi+ (GeV)");
    hiMT31[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiMT31[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiMT32[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dy);
    sprintf (title, "Reduced Transverse Mass of #pi+ (GeV)");
    hiMT32[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiMT32[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "MTpiplus0%s%4.2f", tag1Name,  yv);
    hiMT30.push_back((TH1F*)hiMT32[ii]->Clone());
    hiMT30[ii]->SetName(name); 
    hiMT30[ii]->Scale(scale);
    hiMT30[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin  = hiMT41[ii]->GetBinWidth(1);
    sprintf (title, "Reduced Transverse Mass of #pi- (GeV)");
    hiMT41[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiMT41[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiMT42[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dy);
    sprintf (title, "Reduced Transverse Mass of #pi- (GeV)");
    hiMT42[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiMT42[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "MTpiminus0%s%4.2f", tag1Name,  yv);
    hiMT40.push_back((TH1F*)hiMT42[ii]->Clone());
    hiMT40[ii]->SetName(name); 
    hiMT40[ii]->Scale(scale);
    hiMT40[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin  = hiMT51[ii]->GetBinWidth(1);
    sprintf (title, "Reduced Transverse Mass of K+ (GeV)");
    hiMT51[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiMT51[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiMT52[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dy);
    sprintf (title, "Reduced Transverse Mass of K+ (GeV)");
    hiMT52[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiMT52[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "MTkplus0%s%4.2f", tag1Name,  yv);
    hiMT50.push_back((TH1F*)hiMT52[ii]->Clone());
    hiMT50[ii]->SetName(name); 
    hiMT50[ii]->Scale(scale);
    hiMT50[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin  = hiMT61[ii]->GetBinWidth(1);
    sprintf (title, "Reduced Transverse Mass of K- (GeV)");
    hiMT61[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events/%6.3f GeV", xbin);
    hiMT61[ii]->GetYaxis()->SetTitle(title);
    xbin  = hiMT62[ii]->GetBinWidth(1);
    scale = cross_sec/(((double)(std::max(nevt,1)))*xbin*2.*pi*dy);
    sprintf (title, "Reduced Transverse Mass of K- (GeV)");
    hiMT62[ii]->GetXaxis()->SetTitle(title);
    sprintf (title, "Events (scaled by #frac{1}{p})/%6.3f GeV", xbin);
    hiMT62[ii]->GetYaxis()->SetTitle(title);
    sprintf (name, "MTkminus0%s%4.2f", tag1Name,  yv);
    hiMT60.push_back((TH1F*)hiMT62[ii]->Clone());
    hiMT60[ii]->SetName(name); 
    hiMT60[ii]->Scale(scale);
    hiMT60[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");
  }

  TFile f(fileName.c_str(), "recreate");
  for (unsigned int ii=0; ii<ymin.size(); ii++) {
    hiMT11[ii]->Write(); hiMT10[ii]->Write(); hiMT12[ii]->Write();
    hiMT21[ii]->Write(); hiMT20[ii]->Write(); hiMT22[ii]->Write();
    hiMT31[ii]->Write(); hiMT30[ii]->Write(); hiMT32[ii]->Write();
    hiMT41[ii]->Write(); hiMT40[ii]->Write(); hiMT42[ii]->Write();
    hiMT51[ii]->Write(); hiMT50[ii]->Write(); hiMT52[ii]->Write();
    hiMT61[ii]->Write(); hiMT60[ii]->Write(); hiMT62[ii]->Write();
  }
  epTest.write();
  f.Close();
}

void HistoBNLTest47::initialize() {

  unInitialized = false;
  G4cout << "HistoBNLTest47::initialize invoked" << G4endl;
  fileName = particle + target + generator;
  std::ostringstream tmp;
  char nams[4];
  sprintf (nams, "%4.1f", energy);
  tmp << nams << "GeV";
  if ( jobID > -1 )     tmp << "-" << jobID;
  if ( clusterID > -1 ) tmp << "-" << clusterID;
  tmp << ".root";
  fileName += tmp.str();
  sprintf (tag1Name, "%s%s%s%4.1fGeV", particle.c_str(), target.c_str(),
	   generator.c_str(), energy); 
  sprintf (tag2Name, "%s+%s", particle.c_str(), target.c_str());
  sprintf (tag3Name, "at %4.1f GeV (%s)", energy, generator.c_str());
  G4cout << "HistoBNLTest47::fileName:" << fileName.c_str() 
	 << " Tag1:" << tag1Name << " Tag2: " << tag2Name << " Tag3: " 
	 << tag3Name << G4endl;

  book();
  epTest.initialize(particle,target,energy,generator);
}

void HistoBNLTest47::book() {

  char name[100], title[160];
  hiMT11.clear(); hiMT12.clear(); hiMT21.clear(); hiMT22.clear();
  hiMT31.clear(); hiMT32.clear(); hiMT41.clear(); hiMT42.clear();
  hiMT51.clear(); hiMT52.clear(); hiMT61.clear(); hiMT62.clear();
  for (unsigned int ii=0; ii<rapidities.size()-1; ii++) {
    double yv = 0.5*(rapidities[ii]+rapidities[ii+1]);
    sprintf (name, "MTproton1%s%4.2f", tag1Name,  yv);
    sprintf (title, "%s to p %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT11.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTproton2%s%4.2f", tag1Name,  yv);
    hiMT12.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTneutron1%s%4.2f", tag1Name,  yv);
    sprintf (title, "%s to n %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT21.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTneutron2%s%4.2f", tag1Name,  yv);
    hiMT22.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTpiplus1%s%4.2f", tag1Name,  yv);
    sprintf (title, "%s to #pi+ %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT31.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTpiplus2%s%4.2f", tag1Name,  yv);
    hiMT32.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTpiminus1%s%4.2f", tag1Name,  yv);
    sprintf (title, "%s to #pi- %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT41.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTpiminus2%s%4.2f", tag1Name,  yv);
    hiMT42.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTkplus1%s%4.2f", tag1Name,  yv);
    sprintf (title, "%s to K+ %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT51.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTkplus2%s%4.2f", tag1Name,  yv);
    hiMT52.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTkminus1%s%4.2f", tag1Name,  yv);
    sprintf (title, "%s to K- %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT61.push_back(new TH1F (name, title, 800, 0., 8.));
    sprintf (name, "MTkminus2%s%4.2f", tag1Name,  yv);
    hiMT62.push_back(new TH1F (name, title, 800, 0., 8.));
  }
}
