#include "mkpldda.h"

#include "mkplexpdatamanager.h"
#include "TPostScript.h"

const int mkpldda::mkpl_num_bins = 100;


void mkpldda::PrepareHistograms(mkplexpdatamanager * expdata)
{
  // Initialize histograms
  this->InitializeHistograms();
  
  // Get experimental data
  mkpl_exp = expdata->GetDDA(mkpl_A, mkpl_Z, mkpl_projectile_A, mkpl_projectile_Z,
			     mkpl_target_A, mkpl_target_Z, mkpl_reaction_E, 
			     !this->IsCM()); // False for CM
                                             // True  for LAB

  if (mkpl_exp)
    {
      TString afijo;
      if (this->IsCM())
	{
	  afijo = "_cm_dda_";
	}
      else 
	{
	  afijo = "_dda_";
	}
      for (vector<mkplexpdda>::iterator it = mkpl_exp->begin();
	   it != mkpl_exp->end(); it++)
	{
	  double rang_l = it->GetInfCut();
	  double rang_h = it->GetSupCut();
	  TString tag("mkpl_");
	  TString tagpreeq("mkpl_preeq_");
	  TString tagevap("mkpl_evap_");
	  TString tagfis("mkpl_fis_");
	  TString tagfermi("mkpl_fermi_");
	  TString taginc("mkpl_inc_");
	  tag      += mkpl_name + afijo;
	  tagpreeq += mkpl_name + afijo;
	  tagevap  += mkpl_name + afijo;
	  tagfis   += mkpl_name + afijo;
	  tagfermi += mkpl_name + afijo;
	  taginc   += mkpl_name + afijo;
	  tag      += rang_l;
	  tagpreeq += rang_l;
	  tagevap  += rang_l;
	  tagfis   += rang_l;
	  tagfermi += rang_l;
	  taginc   += rang_l;
	  tag      += rang_h;
	  tagpreeq += rang_h;
	  tagevap  += rang_h;
	  tagfis   += rang_h;
	  tagfermi += rang_h;
	  taginc   += rang_h;
	  tag.ReplaceAll(" ","");
	  tagpreeq.ReplaceAll(" ","");
	  tagevap.ReplaceAll(" ","");
	  tagfis.ReplaceAll(" ","");
	  tagfermi.ReplaceAll(" ","");
	  taginc.ReplaceAll(" ","");

	  
	  mkpl_h.push_back(new TH1F(tag, it->GetData()->GetTitle(),mkpl_num_bins, 
				    0.0, 180.0));
      
	  mkpl_preeq_h.push_back(new TH1F(tagpreeq, it->GetData()->GetTitle(), mkpl_num_bins, 
					  0.0, 180.0)); 
	  
	  mkpl_evap_h.push_back(new TH1F(tagevap, it->GetData()->GetTitle(), mkpl_num_bins, 
					 0.0, 180.0)); 
      
	  mkpl_fis_h.push_back(new TH1F(tagfis, it->GetData()->GetTitle(), mkpl_num_bins, 
					0.0, 180.0)); 
      
	  mkpl_fermi_h.push_back(new TH1F(tagfermi, it->GetData()->GetTitle(), mkpl_num_bins, 
					  0.0, 180.0)); 
      
	  mkpl_inc_h.push_back(new TH1F(taginc, it->GetData()->GetTitle(), mkpl_num_bins, 
					0.0, 180.0)); 
	}
    }
  return;

}


void mkpldda::Normalization(const double xs, const double entries)
{
  mkpl_weight.clear();
  if (!mkpl_h.empty())
    {
      for (unsigned int i = 0; i < mkpl_exp->size(); i++)
	{
	  double binwidth = mkpl_h[i]->GetBinWidth(1)*degtorad;
	  double norma = xs / (entries * 2.0*TMath::Pi() * binwidth );
	  mkpl_weight.push_back(norma);
	}
    }
  return;
}


void mkpldda::Fill(const double E, const double ang, const string& pname)
{
  if (mkpl_exp)
    {
      for (vector<mkplexpdda>::iterator it = mkpl_exp->begin();
	   it != mkpl_exp->end(); it++)
	{
	  if (E >= it->GetInfCut() && E <= it->GetSupCut())
	    {
	      int n = distance(mkpl_exp->begin(),it);
	      double ang_deg = ang * radtodeg;
	      double sinus = TMath::Sin(mkpl_h[n]->GetBinCenter(mkpl_h[n]->FindBin(ang_deg))*degtorad);
	      mkpl_h[n]->Fill(ang_deg,mkpl_weight[n]/sinus);
	      if (pname == mkpl_preeq_model_name)
		{
		  mkpl_preeq_h[n]->Fill(ang_deg, mkpl_weight[n]/sinus);
		}
	      else if (pname == mkpl_evap_model_name)
		{
		  mkpl_evap_h[n]->Fill(ang_deg, mkpl_weight[n]/sinus);
		} 
	      else if (pname == mkpl_fis_model_name)
		{
		  mkpl_fis_h[n]->Fill(ang_deg, mkpl_weight[n]/sinus);
		} 
	      else if (pname == mkpl_fermi_model_name)
		{
		  mkpl_fermi_h[n]->Fill(ang_deg, mkpl_weight[n]/sinus);
		} 
	      else if (pname == mkpl_inc_model_name)
		{
		  mkpl_inc_h[n]->Fill(ang_deg, mkpl_weight[n]/sinus);
		} 
	    }
	}
    }
  return;
}

vector<mkplcomparisonagreement> mkpldda::Plot(bool comp, const int n)
{
  if (!mkpl_h.empty())
    {
      if (!this->IsCM())
	{
	  mkpl_h[n]->GetXaxis()->SetTitle("#theta (deg)");
	  mkpl_h[n]->GetYaxis()->SetTitle("#frac{d#sigma}{d#Theta} (mb/rad)");
	}
      else 
	{
	  mkpl_h[n]->GetXaxis()->SetTitle("#theta (MeV) (CM)");
	  mkpl_h[n]->GetYaxis()->SetTitle("#frac{d#sigma}{d#Theta} (mb/rad) (CM)");
	}
      double maxExp = 0.0;
      double minExp = DBL_MAX;
      int color = 2;
      int style = 20;

      for (int j = 0; j < mkpl_exp->operator[](n).GetData()->
	     GetListOfGraphs()->GetSize(); j++)  // Loop for graphs
	{
	  TGraph * g = (TGraph*)(mkpl_exp->operator[](n).GetData()->GetListOfGraphs()->At(j));
	  double * Y = g->GetY();
	  for (int k = 0; k < g->GetN(); k++)
	    {
	      if (Y[k] > maxExp) maxExp = Y[k];
	      if (Y[k] < minExp && Y[k] > 0.0) minExp = Y[k];
	    }
	  g->SetMarkerColor(color++);
	  g->SetMarkerStyle(style++);
	  g->SetMarkerSize(0.75);
	}
      double maxSim = mkpl_h[n]->GetBinContent(mkpl_h[n]->GetMaximumBin());
      double minSim = this->GetMinimumNonNull(mkpl_h[n]);

      maxExp = TMath::Max(maxExp,maxSim);
      minExp = TMath::Min(minExp,minSim);
      if (minExp <= 0.0) minExp = 0.0001;
      else minExp = 0.9*minExp;

      mkpl_h[n]->SetMinimum(minExp);
      mkpl_h[n]->SetMaximum(1.2*maxExp);
      mkpl_h[n]->Draw();
      if (comp)
	{
	  mkpl_preeq_h[n]->SetLineColor(6);
	  mkpl_preeq_h[n]->Draw("SAME");
	  mkpl_evap_h[n]->SetLineColor(7);
	  mkpl_evap_h[n]->Draw("SAME");
	  mkpl_fis_h[n]->SetLineColor(8);
	  mkpl_fis_h[n]->Draw("SAME");
	  mkpl_fermi_h[n]->SetLineColor(30);
	  mkpl_fermi_h[n]->Draw("SAME");
	  mkpl_inc_h[n]->SetLineColor(9);
	  mkpl_inc_h[n]->Draw("SAME");
	}
      mkpl_h[n]->Draw("SAME");	
      mkpl_exp->operator[](n).GetData()->Draw("P"); 
    }
  return this->Agreement(mkpl_exp->operator[](n).GetData(), mkpl_h[n]);
}


void mkpldda::SaveHistograms()
{
  for (unsigned int i = 0; i < mkpl_exp->size(); i++)
    {
      mkpl_h[i]->Write(Form("simulation_range_%u",i));	  
      mkpl_preeq_h[i]->Write(Form("simpreeq_range_%u",i));	  
      mkpl_evap_h[i]->Write(Form("simevap_range_%u",i));	  
      mkpl_fermi_h[i]->Write(Form("simfermi_range_%u",i));	  
      mkpl_inc_h[i]->Write(Form("siminc_range_%u",i));	  
      mkpl_exp->operator[](i).GetData()->Write(Form("experimental_angle_%u",i));
    }
  return;
}


void mkpldda::SaveEPS(const TString& bn, bool comp)
{
  TString fnamebase(bn);
  fnamebase += Form("-%u-%u-",mkpl_A,mkpl_Z);
  if (this->IsCM()) fnamebase += "-CM-";
  fnamebase += "DDA-r";
  TPostScript epsfile;
  for (unsigned int i = 0; i < mkpl_exp->size(); i++)
    {
      TString fname = fnamebase;
      fname += i;
      fname += ".eps";
      epsfile.Open(fname,113);
      epsfile.Range(20.0,20.0);
      this->Plot(comp,i);
      epsfile.Close();
    }
  return;
}
