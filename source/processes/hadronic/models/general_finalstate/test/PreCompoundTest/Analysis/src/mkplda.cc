#include "mkplda.h"

#include "mkplexpdatamanager.h"
#include "TPostScript.h"

const int mkplda::mkpl_num_bins = 100;


void mkplda::PrepareHistograms(mkplexpdatamanager * expdata)
{
  // Initialize histograms
  this->InitializeHistograms();
  
  // Get experimental data
  mkpl_exp = expdata->GetDA(mkpl_A, mkpl_Z, mkpl_projectile_A, mkpl_projectile_Z,
			    mkpl_target_A, mkpl_target_Z, mkpl_reaction_E, 
			    !this->IsCM()); // False for CM
                                            // True  for LAB

  if (mkpl_exp)
    {
      TString afijo;
      if (this->IsCM())
	{
	  afijo = "_cm_da";
	}
      else 
	{
	  afijo = "_da";
	}
      TString tag = "mkpl_" + mkpl_name + afijo;
      mkpl_h = new TH1F(tag, mkpl_exp->GetTitle(), mkpl_num_bins, 0.0, 180.0);
      
      tag = "mkpl_preeq_" + mkpl_name + afijo;
      mkpl_preeq_h = new TH1F(tag, mkpl_exp->GetTitle(), mkpl_num_bins, 0.0, 180.0); 

      tag = "mkpl_evap_" + mkpl_name + afijo;
      mkpl_evap_h = new TH1F(tag, mkpl_exp->GetTitle(), mkpl_num_bins, 0.0, 180.0); 
      
      tag = "mkpl_fis_" + mkpl_name + afijo;
      mkpl_fis_h = new TH1F(tag, mkpl_exp->GetTitle(), mkpl_num_bins, 0.0, 180.0); 
      
      tag = "mkpl_fermi_" + mkpl_name + afijo;
      mkpl_fermi_h = new TH1F(tag, mkpl_exp->GetTitle(), mkpl_num_bins, 0.0, 180.0); 
      
      tag = "mkpl_inc_" + mkpl_name + afijo;
      mkpl_inc_h = new TH1F(tag, mkpl_exp->GetTitle(), mkpl_num_bins, 0.0, 180.0); 
 
    }
  
  return;

}


void mkplda::Normalization(const double xs, const double entries)
{
  mkpl_weight = 1.0;
  if (mkpl_h)
    {
      double binwidth  = mkpl_h->GetBinWidth(1)*degtorad;
      mkpl_weight = xs / ( entries * 2.0 * TMath::Pi() * binwidth );
    }
  return;
}


void mkplda::Fill(const double E, const double ang, const string& pname)
{
  if (mkpl_h)
    {
      double ang_deg = ang * radtodeg;
      double sinus = TMath::Sin(mkpl_h->GetBinCenter(mkpl_h->FindBin(ang_deg))*degtorad);
      mkpl_h->Fill(ang_deg,mkpl_weight/sinus);
      if (pname == mkpl_preeq_model_name)
	{
	  mkpl_preeq_h->Fill(ang_deg,mkpl_weight/sinus);
	}
      else if (pname == mkpl_evap_model_name)
	{
	  mkpl_evap_h->Fill(ang_deg,mkpl_weight/sinus);
	} 
      else if (pname == mkpl_fermi_model_name)
	{
	  mkpl_fermi_h->Fill(ang_deg,mkpl_weight/sinus);
	} 
      else if (pname == mkpl_fis_model_name)
	{
	  mkpl_fis_h->Fill(ang_deg,mkpl_weight/sinus);
	} 
      else if (pname == mkpl_inc_model_name)
	{
	  mkpl_inc_h->Fill(ang_deg,mkpl_weight/sinus);
	} 
    }
  return;
}

vector<mkplcomparisonagreement> mkplda::Plot(bool comp, const int n)
{
  if (mkpl_exp)
    {

      if (!this->IsCM())
	{
	  mkpl_h->GetXaxis()->SetTitle("#theta (deg)");
	  mkpl_h->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}  (mb/deg)");
	}
      else
	{
	  mkpl_h->GetXaxis()->SetTitle("#theta (deg) (CM)");
	  mkpl_h->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}  (mb/deg) (CM)");
	}
      double maxExp = 0.0;
      double minExp = DBL_MAX;
      int color = 2;
      int style = 20;
      for (int i = 0; i < mkpl_exp->GetListOfGraphs()->GetSize(); i++)
	{
	  TGraph * g = (TGraph*)(mkpl_exp->GetListOfGraphs()->At(i)); 
	  double * Y = g->GetY();
	  for (int j = 0; j < g->GetN(); j++)
	    {
	      if (Y[j] > maxExp) maxExp = Y[j];
	      if (Y[j] < minExp && Y[j] > 0.0) minExp = Y[j];
	    }
	  g->SetMarkerColor(color++);
	  g->SetMarkerStyle(style++);
	  g->SetMarkerSize(0.75);
	}
      double maxSim = mkpl_h->GetBinContent(mkpl_h->GetMaximumBin());
      double minSim = this->GetMinimumNonNull(mkpl_h);

      maxExp = TMath::Max(maxExp,maxSim);
      minExp = TMath::Min(minExp,minSim);
      if (minExp <= 0.0) minExp = 0.0001;
      else minExp = 0.9*minExp;
      
      mkpl_h->SetMinimum(minExp);
      mkpl_h->SetMaximum(1.3*maxExp);
      mkpl_h->Draw();
      if (mkpl_preeq_h && comp) 
	{
	  mkpl_preeq_h->SetLineColor(6);
	  mkpl_preeq_h->Draw("SAME"); 
	  
	  mkpl_evap_h->SetLineColor(7);
	  mkpl_evap_h->Draw("SAME"); 
	  
	  mkpl_fis_h->SetLineColor(8);
	  mkpl_fis_h->Draw("SAME"); 
	  
	  mkpl_fermi_h->SetLineColor(30);
	  mkpl_fermi_h->Draw("SAME"); 
	  
	  mkpl_inc_h->SetLineColor(9);
	  mkpl_inc_h->Draw("SAME"); 
	}
      mkpl_h->Draw("SAME");
      mkpl_exp->Draw("P");
    }
  return this->Agreement(mkpl_exp,mkpl_h);
}

void mkplda::SaveEPS(const TString& bn, bool comp)
{
  TString fname(bn);
  fname += Form("-%u-%u-",mkpl_A,mkpl_Z);
  if (this->IsCM()) fname += "-CM-";
  fname += "DA.eps";
  TPostScript epsfile(fname,113);
  epsfile.Range(20.0,20.0);
  this->Plot(comp);
  epsfile.Close();
  return;
}

