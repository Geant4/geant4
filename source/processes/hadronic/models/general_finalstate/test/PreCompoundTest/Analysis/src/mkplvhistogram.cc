#include "mkplvhistogram.h"

#include "TH1F.h"
#include "TGraph.h"
#include "TMultiGraph.h"

#include "mkplcomparisonagreement.h"

const char * mkplvhistogram::mkpl_preeq_model_name = "G4PreCompoundModel";
const char * mkplvhistogram::mkpl_evap_model_name  = "G4Evaporation";
const char * mkplvhistogram::mkpl_fermi_model_name   = "G4FermiBreakUp";
const char * mkplvhistogram::mkpl_fis_model_name   = "G4CompetitiveFission";
const char * mkplvhistogram::mkpl_inc_model_name   = "G4BinaryCascade";

const double mkplvhistogram::degtorad = TMath::Pi()/180.0;
const double mkplvhistogram::radtodeg = 180.0/TMath::Pi();



vector<mkplcomparisonagreement> mkplvhistogram::Agreement(TMultiGraph * exp, TH1F * sim)
{
  vector<mkplcomparisonagreement> result;
  for (int i = 0; i < exp->GetListOfGraphs()->GetSize(); i++)
    { 
      TGraph * g = (TGraph*)(exp->GetListOfGraphs()->At(i)); 
      int Ns(0);
      double * Sigma_exp = g->GetY();      
      double * Ei = g->GetX();
      double Fmax = 0.0;
      double Fmin = 1.e10;
      double Fbar = 0.0;
      double Favg = 0.0;
      for (int j = 0; j < g->GetN(); j++)
	{
	  // Determine bin for each Ei
	  int bin = sim->FindBin(Ei[j]);
	  double E1 = sim->GetBinCenter(bin);
	  double S1 = sim->GetBinContent(bin);
	  double E2;
	  double S2;
	  if (Ei[j] > E1) 
	    {
	      E2 = sim->GetBinCenter(bin+1);
	      S2 = sim->GetBinContent(bin+1);
	    }
	  else 
	    {
	      E2 = sim->GetBinCenter(bin-1);
	      S2 = sim->GetBinContent(bin-1);
	    }
	  if (S1 > 0.0 && S2 > 0.0)
	    {
	      double LSigma_exp = TMath::Log10(Sigma_exp[j]);
	      double X = TMath::Log10(Ei[j]);
	      E1 = TMath::Log10(E1);
	      E2 = TMath::Log10(E2);
	      S1 = TMath::Log10(S1);
	      S2 = TMath::Log10(S2);
	      // Double logarithmic interpolation
	      double LSigma_sim = ((X-E1)/(E2-E1)) * S2 + ((X-E2)/(E1-E2)) * S1;
	      double Sigma_sim = TMath::Power(10.0,LSigma_sim);
	      double f = Sigma_exp[j]/Sigma_sim;
	      if (f > Fmax) Fmax = f;
	      if (f < Fmin) Fmin = f;
	      double tmp = LSigma_exp - LSigma_sim;
	      Fbar += tmp;
	      Favg += tmp*tmp;
	      Ns++;
	    }
	}
      if (Ns > 0)
	{
	  Fbar = TMath::Power(10.0, Fbar/double(Ns));
	  Favg = TMath::Power(10.0, TMath::Sqrt(Favg)/double(Ns));
	}
      else
	{
	  Fbar = 0.0;
	  Favg = 0.0;
	  Fmax = 0.0;
	  Fmin = 0.0;
	}
      mkplcomparisonagreement dataset;
      dataset.Ns = Ns;
      dataset.Fbar = Fbar;
      dataset.Favg = Favg;
      dataset.Fmax = Fmax;
      dataset.Fmin = Fmin;
      result.push_back(dataset);
    }
  return result;

}


double mkplvhistogram::GetMinimumNonNull(TH1F * h)
{
  double min = h->GetBinContent(h->GetMinimumBin());
  if (min <= 0.0)
    {
      min = DBL_MAX;
      for (Int_t b = 0; b < h->GetNbinsX(); b++)
	{
	  double tmp = h->GetBinContent(b);
	  if (tmp > 0.0 && tmp < min) min = tmp;
	}
    }
  return min;
}
