#include "mkplcomparisonhistograms.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TStyle.h"

void mkplcomparisonhistograms::DeleteHistograms()
{
  // delete all existing DE histograms
  if (mkpl_de) delete mkpl_de;
  if (mkpl_cm_de) delete mkpl_cm_de;

  // delete all existing DA histograms
  if (mkpl_da) delete mkpl_da;
  if (mkpl_cm_da) delete mkpl_cm_da;

  // delete all existing DD histograms
  if (mkpl_dd) delete mkpl_dd;
  if (mkpl_cm_dd) delete mkpl_cm_dd; 

  // delete all existing DD histograms
  if (mkpl_dda) delete mkpl_dda;
  if (mkpl_cm_dda) delete mkpl_cm_dda;
  return;
}



void mkplcomparisonhistograms::InitializeHistograms()
{
  this->DeleteHistograms();
  mkpl_de = 0;
  mkpl_cm_de = 0;
  mkpl_da = 0;
  mkpl_cm_da = 0;
  mkpl_dd = 0;
  mkpl_cm_dd = 0;
  mkpl_dda = 0;
  mkpl_cm_dda = 0;
  return;
}


void mkplcomparisonhistograms::PrepareDEhistograms(mkplexpdatamanager * expdata,
						   const bool components)
{
  // Get experimental data
  mkpl_de = new mkplde(mkpl_A, mkpl_Z,
		       mkpl_projectile_A, mkpl_projectile_Z,
		       mkpl_target_A, mkpl_target_Z,
		       mkpl_reaction_E, mkpl_name);
  mkpl_de->PrepareHistograms(expdata);

  if (!mkpl_de->ThereIsData())
    {
      delete mkpl_de;
      mkpl_de = 0;
    }

  mkpl_cm_de = new mkplde(mkpl_A, mkpl_Z,
			  mkpl_projectile_A, mkpl_projectile_Z,
			  mkpl_target_A, mkpl_target_Z,
			  mkpl_reaction_E,mkpl_name);
  mkpl_cm_de->SetCM();
  mkpl_cm_de->PrepareHistograms(expdata);

  if (!mkpl_cm_de->ThereIsData())
    {
      delete mkpl_cm_de;
      mkpl_cm_de = 0;
    }
  
  return;
}

void mkplcomparisonhistograms::PrepareDAhistograms(mkplexpdatamanager * expdata,
						   const bool components)
{
  // Get experimental data
  mkpl_da = new mkplda(mkpl_A, mkpl_Z,
		       mkpl_projectile_A, mkpl_projectile_Z,
		       mkpl_target_A, mkpl_target_Z,
		       mkpl_reaction_E, mkpl_name);
  mkpl_da->PrepareHistograms(expdata);

  if (!mkpl_da->ThereIsData())
    {
      delete mkpl_da;
      mkpl_da = 0;
    }

  mkpl_cm_da = new mkplda(mkpl_A, mkpl_Z,
			  mkpl_projectile_A, mkpl_projectile_Z,
			  mkpl_target_A, mkpl_target_Z,
			  mkpl_reaction_E, mkpl_name);
  mkpl_cm_da->SetCM();
  mkpl_cm_da->PrepareHistograms(expdata);

  if (mkpl_cm_da)
    {
      delete mkpl_cm_da;
      mkpl_cm_da = 0;
    }
  return;
}


void mkplcomparisonhistograms::PrepareDDhistograms(mkplexpdatamanager * expdata,
						   const bool components)
{ 
  // Get experimental data
  mkpl_dd = new mkpldd(mkpl_A, mkpl_Z,
		       mkpl_projectile_A, mkpl_projectile_Z,
		       mkpl_target_A, mkpl_target_Z,
		       mkpl_reaction_E, mkpl_name);
  mkpl_dd->PrepareHistograms(expdata);

  if (!mkpl_dd->ThereIsData())
    {
      delete mkpl_dd;
      mkpl_dd = 0;
    }


  mkpl_cm_dd = new mkpldd(mkpl_A, mkpl_Z,
			  mkpl_projectile_A, mkpl_projectile_Z,
			  mkpl_target_A, mkpl_target_Z,
			  mkpl_reaction_E, mkpl_name);

  mkpl_cm_dd->SetCM();
  mkpl_cm_dd->PrepareHistograms(expdata);

  if (!mkpl_cm_dd->ThereIsData())
    {
      delete mkpl_cm_dd;
      mkpl_cm_dd = 0;
    }

  return;
}


void mkplcomparisonhistograms::PrepareDDAhistograms(mkplexpdatamanager * expdata,
						    const bool components)
{ 
  // Get experimental data
  mkpl_dda = new mkpldda(mkpl_A, mkpl_Z,
			 mkpl_projectile_A, mkpl_projectile_Z,
			 mkpl_target_A, mkpl_target_Z,
			 mkpl_reaction_E, mkpl_name);
  mkpl_dda->PrepareHistograms(expdata);

  if (!mkpl_dda->ThereIsData())
    {
      delete mkpl_dda;
      mkpl_dda = 0;
    }
  
  
  mkpl_cm_dda = new mkpldda(mkpl_A, mkpl_Z,
			    mkpl_projectile_A, mkpl_projectile_Z,
			    mkpl_target_A, mkpl_target_Z,
			    mkpl_reaction_E, mkpl_name);
  mkpl_cm_dda->SetCM();
  mkpl_cm_dda->PrepareHistograms(expdata);

  if (!mkpl_cm_dda->ThereIsData())
    {
      delete mkpl_cm_dda;
      mkpl_cm_dda = 0;
    }
  return;
}




void mkplcomparisonhistograms::DENormalization(const double xs, const double entries)
{
  if (mkpl_de)
    {
      mkpl_de->Normalization(xs,entries);
    }
  if (mkpl_cm_de)
    {
      mkpl_cm_de->Normalization(xs,entries);
    }
  return;
}

void mkplcomparisonhistograms::DANormalization(const double xs, const double entries)
{
  if (mkpl_da)
    {
      mkpl_da->Normalization(xs,entries);
    }
  if (mkpl_cm_da)
    {
      mkpl_cm_da->Normalization(xs,entries);
    }
  return;
}


void mkplcomparisonhistograms::DDNormalization(const double xs, const double entries)
{
  if (mkpl_dd)
    {
      mkpl_dd->Normalization(xs,entries);
    }

  if (mkpl_cm_dd)
    {
      mkpl_cm_dd->Normalization(xs,entries);
    }
  return;
}

void mkplcomparisonhistograms::DDANormalization(const double xs, const double entries)
{
  if (mkpl_dda)
    {
      mkpl_dda->Normalization(xs,entries);
    }

  if (mkpl_cm_dda)
    {
      mkpl_cm_dda->Normalization(xs,entries);
    }
  return;
}


void mkplcomparisonhistograms::FillDE(const double E, const double ang, const string& pname)
{
  if (mkpl_de)
    {
      mkpl_de->Fill(E,ang,pname);
    }
  return;
}

void mkplcomparisonhistograms::FillDEcm(const double E, const double ang, const string& pname)
{
  if (mkpl_cm_de)
    {
      mkpl_cm_de->Fill(E,ang,pname);
    }
  return;
}

void mkplcomparisonhistograms::FillDA(const double E, const double ang, const string& pname)
{
  if (mkpl_da)
    {
      mkpl_da->Fill(E,ang,pname);
    }
  return;
}

void mkplcomparisonhistograms::FillDAcm(const double E, const double ang, const string& pname)
{
  if (mkpl_cm_da)
    {
      mkpl_cm_da->Fill(E,ang,pname);
    }
  return;
}


void mkplcomparisonhistograms::FillDD(const double E, const double ang, const string& pname)
{
  if (mkpl_dd)
    {
      mkpl_dd->Fill(E,ang,pname);
    }
  return;
}

void mkplcomparisonhistograms::FillDDcm(const double E, const double ang, const string& pname)
{
  if (mkpl_cm_dd)
    {
      mkpl_cm_dd->Fill(E,ang,pname);
    }
  return;
}

void mkplcomparisonhistograms::FillDDA(const double E, const double ang, const string& pname)
{
  if (mkpl_dda)
    {
      mkpl_dda->Fill(E,ang,pname);
    }
  return;
}

void mkplcomparisonhistograms::FillDDAcm(const double E, const double ang, const string& pname)
{
  if (mkpl_cm_dda)
    {
      mkpl_cm_dda->Fill(E,ang,pname);
    }
  return;
}


void mkplcomparisonhistograms::SaveHistograms(TFile * f)
{
  f->cd();
  TDirectory * rootdir = f->mkdir(mkpl_name);
  rootdir->cd();
  if (mkpl_de || mkpl_cm_de)
    {
      TDirectory * DEdir = rootdir->mkdir("DE","Differential Energy Cross Section");
      DEdir->cd();
      if (mkpl_de)
	{
	  TDirectory * DELABdir = DEdir->mkdir("LAB","LAB reference frame");
	  DELABdir->cd();
	  mkpl_de->SaveHistograms();
	}
      if (mkpl_cm_de)
	{
	  TDirectory * DECMdir = DEdir->mkdir("CM","CM reference frame");
	  DECMdir->cd();
	  mkpl_cm_de->SaveHistograms();
	}
    }
  if (mkpl_da || mkpl_cm_da || mkpl_dda || mkpl_cm_dda)
    {
      TDirectory * DAdir = rootdir->mkdir("DA","Differential Angular Cross Section");
      DAdir->cd();
      if (mkpl_da || mkpl_dda)
	{
	  TDirectory * DALABdir = DAdir->mkdir("LAB","LAB reference frame");
	  DALABdir->cd();
	  if (mkpl_da)
	    {
	      mkpl_da->SaveHistograms();
	    }
	  if (mkpl_dda)
	    {
	      mkpl_dda->SaveHistograms();
	    }
	}
      if (mkpl_cm_da || mkpl_cm_dda)
	{
	  TDirectory * DACMdir = DAdir->mkdir("CM","CM reference frame");
	  DACMdir->cd();      
	  if (mkpl_cm_da)
	    {
	      mkpl_cm_da->SaveHistograms();
	    }
	  if (mkpl_cm_dda)
	    {
	      mkpl_cm_dda->SaveHistograms();
	    }
	}
    }
  if (mkpl_dd || mkpl_cm_dd)
    {
      TDirectory * DDdir = rootdir->mkdir("DD","Double Differential Energy-Angle Cross Section");
      DDdir->cd();
      if (mkpl_dd)
	{
	  TDirectory * DDLABdir = DDdir->mkdir("LAB","LAB reference frame");
	  DDLABdir->cd();
	  mkpl_dd->SaveHistograms();
	}
      if (mkpl_cm_dd)
	{
	  TDirectory * DDCMdir = DDdir->mkdir("CM","CM reference frame");
	  DDCMdir->cd();
	  mkpl_cm_dd->SaveHistograms();
	}
    }
  f->cd();
  return;
}

void mkplcomparisonhistograms::SaveEPS(const TString& bn, bool comp)
{
  if (mkpl_de)
    {
      mkpl_de->SaveEPS(bn,comp);
    }
  if (mkpl_cm_de)
    {
      mkpl_cm_de->SaveEPS(bn,comp);
    }

  if (mkpl_da)
    {
      mkpl_da->SaveEPS(bn,comp);
    }
  if (mkpl_cm_da)
    {
      mkpl_cm_da->SaveEPS(bn,comp);
    }


  if (mkpl_dda)
    {
      mkpl_dda->SaveEPS(bn,comp);
    }
  if (mkpl_cm_dda)
    {
      mkpl_cm_dda->SaveEPS(bn,comp);
    }

  if (mkpl_dd)
    {
      mkpl_dd->SaveEPS(bn,comp);
    }
  if (mkpl_cm_dd)
    {
      mkpl_cm_dd->SaveEPS(bn,comp);
    }

  return;
}



vector<mkplcomparisonagreement> mkplcomparisonhistograms::PlotDE(bool comp)
{
  if (mkpl_de)
    {
      return mkpl_de->Plot(comp);
    }
  return vector<mkplcomparisonagreement>(0);
}


vector<mkplcomparisonagreement> mkplcomparisonhistograms::PlotDEcm(bool comp)
{
  if (mkpl_cm_de)
    {
      return mkpl_cm_de->Plot(comp);
    }
  return vector<mkplcomparisonagreement>(0);
}


vector<mkplcomparisonagreement> mkplcomparisonhistograms::PlotDA(bool comp)
{
  if (mkpl_da)
    {
      return mkpl_da->Plot(comp);
    }
  return vector<mkplcomparisonagreement>(0);
}

vector<mkplcomparisonagreement> mkplcomparisonhistograms::PlotDAcm(bool comp)
{
  if (mkpl_cm_da)
    {
      return mkpl_cm_da->Plot(comp);
    }
  return vector<mkplcomparisonagreement>(0);
}

vector<mkplcomparisonagreement> mkplcomparisonhistograms::PlotDD(const int n, bool comp)
{
  if (mkpl_dd && n > -1 && n <= mkpl_dd->Size())
    {
      return mkpl_dd->Plot(comp,n);
    }
  return vector<mkplcomparisonagreement>(0);
}

vector<mkplcomparisonagreement> mkplcomparisonhistograms::PlotDDcm(const int n, bool comp)
{
  if (mkpl_cm_dd && n > -1 && n <= mkpl_cm_dd->Size())
    {
      return mkpl_cm_dd->Plot(comp,n);
    }
  return vector<mkplcomparisonagreement>(0);
}

vector<mkplcomparisonagreement> mkplcomparisonhistograms::PlotDDA(const int n, bool comp)
{
  if (mkpl_dda && n > -1 && n < mkpl_dda->Size())
    {
      return mkpl_dda->Plot(comp,n);
    }
  return vector<mkplcomparisonagreement>(0);
}

vector<mkplcomparisonagreement> mkplcomparisonhistograms::PlotDDAcm(const int n, bool comp)
{
  if (mkpl_cm_dda && n > -1 && n < mkpl_cm_dda->Size())
    {
      mkpl_cm_dda->Plot(comp,n);
    }
  return vector<mkplcomparisonagreement>(0);
}


