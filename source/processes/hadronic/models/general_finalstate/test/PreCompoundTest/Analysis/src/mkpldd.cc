#include "mkpldd.h"

#include "mkplexpdatamanager.h"
#include "TLatex.h"
#include "TPostScript.h"

const int mkpldd::mkpl_num_bins = 100;
//const double mkpldd::mkpl_angular_window = 2.5 * mkplvhistogram::degtorad;
const double mkpldd::mkpl_angular_window = 0.043; // rad
void mkpldd::PrepareHistograms(mkplexpdatamanager * expdata)
{
  // Initialize histograms
  this->InitializeHistograms();
  
  // Get experimental data
  mkpl_exp = expdata->GetDD(mkpl_A, mkpl_Z, mkpl_projectile_A, mkpl_projectile_Z,
			    mkpl_target_A, mkpl_target_Z, mkpl_reaction_E, 
			    !this->IsCM()); // False for CM
                                            // True  for LAB

  if (mkpl_exp)
    {
      TString afijo;
      if (this->IsCM())
	{
	  afijo = "_cm_dd_";
	}
      else 
	{
	  afijo = "_dd_";
	}
      for (vector<pair<double,TMultiGraph*> >::iterator it = mkpl_exp->begin();
	   it != mkpl_exp->end(); it++)
	{
	  double ang = it->first;
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
	  tag      += ang;
	  tagpreeq += ang;
	  tagevap  += ang;
	  tagfis   += ang;
	  tagfermi += ang;
	  taginc   += ang;
	  tag.ReplaceAll(" ","");
	  tagpreeq.ReplaceAll(" ","");
	  tagevap.ReplaceAll(" ","");
	  tagfis.ReplaceAll(" ","");
	  tagfermi.ReplaceAll(" ","");
	  taginc.ReplaceAll(" ","");

	  
	  mkpl_h.push_back(new TH1F(tag, it->second->GetTitle(),mkpl_num_bins, 
				    0.0, 1.1*mkpl_reaction_E));
      
	  mkpl_preeq_h.push_back(new TH1F(tagpreeq, it->second->GetTitle(), mkpl_num_bins, 
					  0.0, 1.1*mkpl_reaction_E)); 
	  
	  mkpl_evap_h.push_back(new TH1F(tagevap, it->second->GetTitle(), mkpl_num_bins, 
					 0.0, 1.1*mkpl_reaction_E)); 
      
	  mkpl_fis_h.push_back(new TH1F(tagfis, it->second->GetTitle(), mkpl_num_bins, 
					0.0, 1.1*mkpl_reaction_E)); 
      
	  mkpl_fermi_h.push_back(new TH1F(tagfermi, it->second->GetTitle(), mkpl_num_bins, 
					  0.0, 1.1*mkpl_reaction_E)); 
	  
	  mkpl_inc_h.push_back(new TH1F(taginc, it->second->GetTitle(), mkpl_num_bins, 
					0.0, 1.1*mkpl_reaction_E)); 
	}
    }
  return;

}


void mkpldd::Normalization(const double xs, const double entries)
{
  mkpl_weight.clear();
  if (!mkpl_h.empty())
    {
      for (vector<pair<double,TMultiGraph*> >::iterator it = mkpl_exp->begin();
	   it != mkpl_exp->end(); it++)
	{
	  double ang = it->first; // angle is in radians
	  double Thmin = TMath::Max(ang - mkpl_angular_window, 0.0);
	  double Thmax = TMath::Min(ang + mkpl_angular_window, TMath::Pi());
	  if (ang == 0.0)
	    {
	      Thmin = 0.0;
	      Thmax = 2.0 * mkpl_angular_window;
	    }
	  double norma = xs / (entries * mkpl_h[distance(mkpl_exp->begin(),it)]->GetBinWidth(1) *
			       2.0*TMath::Pi() * (TMath::Cos(Thmin)-TMath::Cos(Thmax)));
	  mkpl_weight.push_back(norma);
	}
    }
  return;
}


void mkpldd::Fill(const double E, const double ang, const string& pname)
{
  if (mkpl_exp)
    {
      for (vector<pair<double,TMultiGraph*> >::iterator it = mkpl_exp->begin();
	   it != mkpl_exp->end(); it++)
	{
	  if (TMath::Abs(ang-it->first) <= mkpl_angular_window)
	    {
	      int n = distance(mkpl_exp->begin(),it);
	      mkpl_h[n]->Fill(E,mkpl_weight[n]);
	      if (pname == mkpl_preeq_model_name)
		{
		  mkpl_preeq_h[n]->Fill(E, mkpl_weight[n]);
		}
	      else if (pname == mkpl_evap_model_name)
		{
		  mkpl_evap_h[n]->Fill(E, mkpl_weight[n]);
		} 
	      else if (pname == mkpl_fermi_model_name)
		{
		  mkpl_fermi_h[n]->Fill(E, mkpl_weight[n]);
		} 
	      else if (pname == mkpl_fis_model_name)
		{
		  mkpl_fis_h[n]->Fill(E, mkpl_weight[n]);
		} 
	      else if (pname == mkpl_inc_model_name)
		{
		  mkpl_inc_h[n]->Fill(E, mkpl_weight[n]);
		} 
	    }
	}
    }
      return;
}

vector<mkplcomparisonagreement> mkpldd::Plot(bool comp, const int n)
{
  if (!mkpl_h.empty())
    {
      if (n == int(mkpl_h.size())) 
	{
	  this->PlotAll();
	  return vector<mkplcomparisonagreement>(0);
	}
      if (!this->IsCM())
	{
	  mkpl_h[n]->GetXaxis()->SetTitle("T (MeV)");
	  mkpl_h[n]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#ThetadT} (mb/MeV/rad)");
	}
      else 
	{
	  mkpl_h[n]->GetXaxis()->SetTitle("T (MeV) (CM)");
	  mkpl_h[n]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#ThetadT} (mb/MeV/rad) (CM)");
	}
      double maxExp = 0.0;
      double minExp = DBL_MAX;
      int color = 2;
      int style = 20;

      for (int j = 0; j < mkpl_exp->operator[](n).second->
	     GetListOfGraphs()->GetSize(); j++)  // Loop for graphs
	{
	  TGraph * g = (TGraph*)(mkpl_exp->operator[](n).second->GetListOfGraphs()->At(j));
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

      minExp = TMath::Min(minSim,minExp);
      maxExp = TMath::Max(maxSim,maxExp);
      if (minExp <= 0.0) minExp = 0.001;
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
      mkpl_exp->operator[](n).second->Draw("P");
    }
  return this->Agreement(mkpl_exp->operator[](n).second, mkpl_h[n]);
}


void mkpldd::PlotAll()
{
  // Copy experimental data and simulated histograms
  this->ClearTemp();
  int numangles = mkpl_exp->size();
  double scale = TMath::Power(10.0,numangles-1);
  double maxexp(0.0);
  double maxsim(0.0);
  double minexp(DBL_MAX);
  double minsim(DBL_MAX);
  vector<TString*> text;
  // Copy and scale experimental data
  for (int i = 0; i < numangles; i++) // loop over angles 
    {
      int color = 2;
      int style = 20;
      double angle = mkpl_exp->operator[](i).first;
      text.push_back(new TString(Form("#Theta = %.0f^{o} (#times 10^{%u})",angle*radtodeg,numangles-1-i)));
      TMultiGraph * mg = new TMultiGraph();
      for (int j = 0; j < mkpl_exp->operator[](i).second->GetListOfGraphs()->GetSize(); j++) // loop over
	{                                                                                    // data sets
	  TGraphErrors * g = (TGraphErrors*)(mkpl_exp->operator[](i).second->GetListOfGraphs()->At(j));
	  TGraphErrors * ng = new TGraphErrors(g->GetN(),g->GetX(),g->GetY(),g->GetEX(),g->GetEY());
	  double * sy = ng->GetY();
	  for (int k = 0; k < ng->GetN(); k++)
	    {
	      sy[k] *= scale;
	      if (sy[k] > maxexp) maxexp = sy[k];
	      if (sy[k] < minexp && sy[k] > 0.0) minexp = sy[k];
	    }
	  ng->SetMarkerStyle(style++);
	  ng->SetMarkerColor(color++);
	  ng->SetMarkerSize(0.50);
	  mg->Add(ng);
	}
      exp.push_back(make_pair(angle, mg));
      scale /= 10.0;
    }
  scale = TMath::Power(10.0,numangles-1);;
  // Copy and scale simulated histograms
  for (unsigned int i = 0; i< mkpl_h.size(); i++)
    {
      TH1F * nh = (TH1F*)(mkpl_h[i]->Clone());
      nh->Scale(scale);
      maxsim = TMath::Max(maxsim,nh->GetBinContent(nh->GetMaximumBin()));
      minsim = TMath::Max(minsim,this->GetMinimumNonNull(nh));
      scale /= 10.0;
      sim.push_back(nh);
    }
  
  // Plot the results
  if (!this->IsCM())
    {
      sim[0]->GetXaxis()->SetTitle("T (MeV)");
      sim[0]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#ThetadT} (mb/MeV/rad)");
    }
  else 
    {
      sim[0]->GetXaxis()->SetTitle("T (MeV) (CM)");
      sim[0]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#ThetadT} (mb/MeV/rad) (CM)");
    }

  TLatex label;
  label.SetTextAlign(11);
  label.SetTextSize(0.0170279);
  label.SetTextColor(4);
  minexp = TMath::Min(minexp,minsim);
  if (minexp <= 0.0) minexp = 0.0001;
  else minexp = 0.9*minexp;
  TString tit(sim[0]->GetTitle());
  tit.Resize(tit.Index("MeV",3,0,TString::kExact)+3);
  sim[0]->SetTitle(tit);
  sim[0]->SetMinimum(minexp);
  sim[0]->SetMaximum(1.5*TMath::Max(maxsim,maxexp));
  sim[0]->Draw();
  double x = 0.2 * ( sim[0]->GetXaxis()->GetXmax() - sim[0]->GetXaxis()->GetXmin() );
  double y = TMath::Log10(1.2 * sim[0]->GetBinContent(sim[0]->FindBin(x)));
  label.DrawLatex(x,y,*text[0]);
  for (unsigned int i = 1; i < sim.size(); i++)
    {
      sim[i]->Draw("SAME");
      x *= 1.03;
      y = TMath::Log10(1.4 * sim[i]->GetBinContent(sim[i]->FindBin(x)));
      label.DrawLatex(x,y,*text[i]);
    }
  for (unsigned int i = 0; i < exp.size(); i++)
    {
      exp[i].second->Draw("P");
    }
  
  for_each(text.begin(), text.end(), DeleteHistogramVector());
  return;
}


void mkpldd::SaveHistograms()
{
  for (unsigned int i = 0; i < mkpl_exp->size(); i++)
    {
      mkpl_h[i]->Write(Form("simulation_angle_%u",i));	  
      mkpl_preeq_h[i]->Write(Form("simpreeq_angle_%u",i));	  
      mkpl_evap_h[i]->Write(Form("simevap_angle_%u",i));	  
      mkpl_fermi_h[i]->Write(Form("simfermi_angle_%u",i));	  
      mkpl_inc_h[i]->Write(Form("siminc_angle_%u",i));	  
      mkpl_exp->operator[](i).second->Write(Form("experimental_angle_%u",i));
    }
  return;
}

void mkpldd::SaveEPS(const TString& bn, bool comp)
{
  TString fnamebase(bn);
  fnamebase += Form("-%u-%u-",mkpl_A,mkpl_Z);
  if (this->IsCM()) fnamebase += "-CM-";
  TString fnameall = fnamebase;
  fnamebase += "DD-ang";
  TPostScript epsfile;
  // First plot all angles together
  fnameall += "DD.eps";
  epsfile.Open(fnameall,113);
  epsfile.Range(20.0,20.0);
  this->PlotAll();
  epsfile.Close();
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
