#include "mkpltesthistograms.h"

#include <algorithm>
#include <iostream>

#include "TROOT.h"
#include "TStyle.h"


const int mkpltesthistograms::mkpl_num_bins_test = 100;

void mkpltesthistograms::DeleteHistograms()
{
  // delete all existing histograms
  if (mkpl_baryonconservation) 
    {
      delete mkpl_baryonconservation;
      mkpl_baryonconservation = 0;
    }
  if (mkpl_chargeconservation)
    {
      delete mkpl_chargeconservation;
      mkpl_chargeconservation = 0;
    }
  if (mkpl_energyconservation)
    {
      delete mkpl_energyconservation;
      mkpl_energyconservation = 0;
    }
  if (mkpl_momentumconservation)
    {
      delete mkpl_momentumconservation;
      mkpl_momentumconservation = 0;
    }
  if (mkpl_pxconservation)
    {
      delete mkpl_pxconservation;
      mkpl_pxconservation = 0;
    }
  if (mkpl_pyconservation)
    {
      delete mkpl_pyconservation;
      mkpl_pyconservation = 0;
    }
  if (mkpl_pzconservation)
    {
      delete mkpl_pzconservation;
      mkpl_pzconservation = 0;
    }
  if (mkpl_theta)
    {
      delete mkpl_theta;
      mkpl_theta = 0;
    }
  if (mkpl_thetaprecom)
    {
      delete mkpl_thetaprecom;
      mkpl_thetaprecom = 0;
    }
  if (mkpl_thetaevap)
    {
      delete mkpl_thetaevap;
      mkpl_thetaevap = 0;
    }
  if (mkpl_phi)
    {
      delete mkpl_phi;
      mkpl_phi = 0;
      mkpl_testPhi = false;
      mkpl_direction = 'R';
    }
  if (mkpl_phiprecom)
    {
      delete mkpl_phiprecom;
      mkpl_phiprecom = 0;
    }
  if (mkpl_phievap)
    {
      delete mkpl_phievap;
      mkpl_phievap = 0;
    }
  if (mkpl_phinucleon)
    {
      delete mkpl_phinucleon;
      mkpl_phinucleon = 0;
    }
  if (mkpl_typefragments) 
    {
      delete mkpl_typefragments;
      mkpl_typefragments = 0;
    }
  return;
}


void mkpltesthistograms::PrepareHistograms(const TVector3& IncidentDirection,
					   const set<string>& partt,
					   const set<string>& proct)
{
  this->DeleteHistograms();
  // Histogram for Baryonic number conservation
  mkpl_baryonconservation = new TH1F("mkpl_baryonconservation",
				     "Baryon number conservation",
				     mkpl_num_bins_test, -10.0, 10.0);
  
  // Histogram for Charge conservation
  mkpl_chargeconservation = new TH1F("mkpl_chargeconservation",
				     "Charge conservation test",
				     mkpl_num_bins_test,-1.0,1.0);
  
  // Histogram for Energy Conservation
  mkpl_energyconservation = new TH1F("mkpl_energyconservation",
				     "Energy conservation test",
				     mkpl_num_bins_test,-1.0,1.0);
  
  // Histograms for Momentum Conservation
  mkpl_momentumconservation = new TH1F("mkpl_momentumconservation",
				       "Momentum conservation test",
				       mkpl_num_bins_test, -1.0,1.0);

  mkpl_pxconservation = new TH1F("mkpl_pxconservation",
				 "P_{x} conservation test",
				 mkpl_num_bins_test, -1.0, 1.0);

  mkpl_pyconservation = new TH1F("mkpl_pyconservation",
				 "P_{y} conservation test",
				 mkpl_num_bins_test, -1.0, 1.0);

  mkpl_pzconservation = new TH1F("mkpl_pzconservation",
				 "P_{z} conservation test",
				 mkpl_num_bins_test, -1.0, 1.0);
    
    
  // Histograms for Angular tests
  mkpl_theta = new TH1F("mkpl_theta",
			"Angle w.r.t. the projectile direction",
			mkpl_num_bins_test, 0.0, TMath::Pi());

  mkpl_thetaprecom = new TH1F("mkpl_thetaprecom",
			      "Angle w.r.t. the projectile direction (preequilibrium)",
			      mkpl_num_bins_test, 0.0, TMath::Pi());

  mkpl_thetaevap = new TH1F("mkpl_thetaevap",
			    "Angle w.r.t. the projectile direction (Evaporation)",
			    mkpl_num_bins_test, 0.0, TMath::Pi());
  
  this->TestPhi(IncidentDirection);
  if (mkpl_testPhi)
    {
      mkpl_phi = new TH1F("mkpl_phi",
			  "Azimutal angle",
			  mkpl_num_bins_test, -TMath::Pi(), TMath::Pi());
      
      mkpl_phiprecom = new TH1F("mkpl_phiprecom",
				"Azimutal angle for Precompound fragments",
				mkpl_num_bins_test, -TMath::Pi(), TMath::Pi());
	
      mkpl_phievap = new  TH1F("mkpl_phiev",
			       "Azimutal angle for Evaporation fragments",
			       mkpl_num_bins_test, -TMath::Pi(), TMath::Pi());
	
      mkpl_phinucleon = new TH1F("mkpl_phinucleon",
				 "Azimutal angle for Nucleons",
				 mkpl_num_bins_test, -TMath::Pi(), TMath::Pi());
    }

  mkpl_particletypes.clear();
  mkpl_processtypes.clear();
  // Copy the contents
  mkpl_particletypes.reserve(partt.size());
  mkpl_particletypes.insert(mkpl_particletypes.begin(),partt.begin(),partt.end()); 
  mkpl_processtypes.reserve(proct.size());
  mkpl_processtypes.insert(mkpl_processtypes.begin(),proct.begin(),proct.end());
  
  // Histograms for kind of particles
  mkpl_typefragments = new TH1F("mkpl_typesfragments",
				"Types of fragments",
				mkpl_particletypes.size(),0.0,
				double(mkpl_particletypes.size()));
  mkpl_typefragments->SetBit(TH1::kCanRebin);
    
  return;
}


void mkpltesthistograms::TestPhi(const TVector3& IncidentDirection)
{
  // Determine whether Phi angle should be tested or not 
  const TVector3 dirX(1.0,0.0,0.0);
  const TVector3 dirY(0.0,1.0,0.0);
  const TVector3 dirZ(0.0,0.0,1.0);
  mkpl_testPhi = false;
  mkpl_direction = 'R';
  if (IncidentDirection == dirX)
    {
      mkpl_testPhi = true;
      mkpl_direction = 'X';
    }
  else if (IncidentDirection == dirY)
    {
      mkpl_testPhi = true;
      mkpl_direction = 'Y';
    }
  else if (IncidentDirection == dirZ)
    {
      mkpl_testPhi = true;
      mkpl_direction = 'Z';
    }
  if (mkpl_testPhi)
    std::cout << "Shooting was done in " << mkpl_direction 
              << " direction: Testing for Azimuthal distribution\n";
  else 
    std::cout << "Shooting was done in a random direction: "
              << "Not testing for Azimuthal distribution\n";
  return;
}

void mkpltesthistograms::FillConservationHistograms(const int deltaA, 
						    const int deltaZ, 
						    const double deltaE, 
						    const double deltaP, 
						    const TVector3 & MomentumTest)
{
  // Test conservation histograms
  if (mkpl_baryonconservation)
    {
      mkpl_baryonconservation->Fill(deltaA);
      mkpl_chargeconservation->Fill(deltaZ);
      mkpl_energyconservation->Fill(deltaE);
      mkpl_momentumconservation->Fill(deltaP);
      mkpl_pxconservation->Fill(MomentumTest.Px());
      mkpl_pyconservation->Fill(MomentumTest.Py());
      mkpl_pzconservation->Fill(MomentumTest.Pz());
    }

  return;
}


void mkpltesthistograms::FillHistograms(const double Theta, 
					const TVector3 & fragP, 
					const string & procname,
					const string & fname)
{
  // Angle tests (only if required)
  if (mkpl_theta) mkpl_theta->Fill(Theta);
  if (mkpl_thetaprecom && procname == "G4PreCompoundModel") mkpl_thetaprecom->Fill(Theta);
  if (mkpl_thetaevap && procname == "G4Evaporation") mkpl_thetaevap->Fill(Theta);
  if (mkpl_testPhi && fragP != TVector3(0.0,0.0,0.0))
    {
      TVector3 fragPr(0.0,0.0,0.0);
      if (mkpl_direction == 'X')
	{
	  fragPr.SetX(fragP.Py());
	  fragPr.SetY(fragP.Pz());
	  fragPr.SetZ(fragP.Px());
	}
      else if (mkpl_direction == 'Y')
	{
	  fragPr.SetX(fragP.Pz());
	  fragPr.SetY(fragP.Px());
	  fragPr.SetZ(fragP.Py());
	}
      else if (mkpl_direction == 'Z')
	{
	  fragPr = fragP;
	}
      double Phi = fragPr.Phi();
      mkpl_phi->Fill(Phi);
      if (procname == "G4Evaporation")
	mkpl_phievap->Fill(Phi);
      else if (procname == "G4PreCompoundModel")
	mkpl_phiprecom->Fill(Phi);
      if (fname == "n" || fname == "p")
	mkpl_phinucleon->Fill(Phi);
      
    }
	
  // Particle type histogram 
  if (mkpl_typefragments)
    {
      vector<string>::iterator ipos = find(mkpl_particletypes.begin(),mkpl_particletypes.end(),fname);
      if (ipos != mkpl_particletypes.end())
	mkpl_typefragments->Fill(ipos->c_str(),1);
    }
  return;
}


void mkpltesthistograms::Draw(TApplication * app)
{
  TCanvas * hCanvas = this->RenewCanvas(0);
  gROOT->SetStyle("Plain");

  // Plot conservation histograms (only if required)
  if (mkpl_baryonconservation)
    {
      gROOT->GetStyle("Plain")->SetOptStat(Int_t(1111111));
      
      // Draw Baryonic Number Conservation Histogram
      mkpl_baryonconservation->GetXaxis()->SetTitle("#DeltaA");
      mkpl_baryonconservation->Draw();
      hCanvas->Update();
      app->Run(kTRUE);
      
      // Draw Charge Conservation Histogram
      hCanvas = RenewCanvas(hCanvas);
      mkpl_chargeconservation->GetXaxis()->SetTitle("#DeltaZ");
      mkpl_chargeconservation->Draw();
      hCanvas->Update();
      app->Run(kTRUE);

      // Draw Energy Conservation Histogram
      hCanvas = RenewCanvas(hCanvas);
      mkpl_energyconservation->GetXaxis()->SetTitle("#Delta E (MeV)");
      mkpl_energyconservation->Draw();
      hCanvas->Update();
      app->Run(kTRUE);

      // Draw Momentum Conservation Histograms
      hCanvas = RenewCanvas(hCanvas);
      mkpl_momentumconservation->GetXaxis()->SetTitle("#Delta P (MeV)");
      mkpl_momentumconservation->Draw();
      hCanvas->Update();
      app->Run(kTRUE);

      hCanvas = RenewCanvas(hCanvas);
      mkpl_pxconservation->GetXaxis()->SetTitle("#Delta P_{x} (MeV)");
      mkpl_pxconservation->Draw();
      hCanvas->Update();
      app->Run(kTRUE);

      hCanvas = RenewCanvas(hCanvas);
      mkpl_pyconservation->GetXaxis()->SetTitle("#Delta P_{y} (MeV)");
      mkpl_pyconservation->Draw();
      hCanvas->Update();
      app->Run(kTRUE);
	
      hCanvas = RenewCanvas(hCanvas);
      mkpl_pzconservation->GetXaxis()->SetTitle("#Delta P_{z} (MeV)");
      mkpl_pzconservation->Draw();
      hCanvas->Update();
      app->Run(kTRUE);

      // Draw Theta Angle Histogram
      hCanvas = RenewCanvas(hCanvas);
      mkpl_theta->GetXaxis()->SetTitle("#theta (rad)");
      mkpl_theta->Draw();
      hCanvas->Update();
      app->Run(kTRUE);

      if (mkpl_testPhi)	
	{
	  hCanvas = RenewCanvas(hCanvas);
	  mkpl_phi->GetXaxis()->SetTitle("#phi (rad)");
	  mkpl_phi->Draw();
	  hCanvas->Update();
	  app->Run(kTRUE);
  
	  hCanvas = RenewCanvas(hCanvas);
	  mkpl_phiprecom->GetXaxis()->SetTitle("#phi (rad)");
	  mkpl_phiprecom->Draw();
	  hCanvas->Update();
	  app->Run(kTRUE);
	    
	  hCanvas = RenewCanvas(hCanvas);
	  mkpl_phievap->GetXaxis()->SetTitle("#phi (rad)");
	  mkpl_phievap->Draw();
	  hCanvas->Update();
	  app->Run(kTRUE);
	    
	  hCanvas = RenewCanvas(hCanvas);
	  mkpl_phinucleon->GetXaxis()->SetTitle("#phi (rad)");
	  mkpl_phinucleon->Draw();
	  hCanvas->Update();
	  app->Run(kTRUE);
	    
	}
	
      // Draw Particle Type Histograms
      hCanvas = RenewCanvas(hCanvas);
      mkpl_typefragments->GetYaxis()->SetTitle("N. of Fragments");
      mkpl_typefragments->LabelsDeflate("X");
      mkpl_typefragments->LabelsOption(">v","X"); 
      mkpl_typefragments->Draw();
      hCanvas->Update();
      app->Run(kTRUE);
    }

  return;
}
