#include "mkplexpdatamanager.h"


void mkplexpdatamanager::Initialize(mkplfilemanager * fm)
{
  mkpl_fm = fm;
  // Experimental Branches
  if (mkpl_DEdata) delete mkpl_DEdata;
  mkpl_DEdata = new expdata_de();
  if (mkpl_DAdata) delete mkpl_DAdata;
  mkpl_DAdata = new expdata_da();
  if (mkpl_DDdata) delete mkpl_DDdata;
  mkpl_DDdata = new expdata_dd();
  if (mkpl_DDAdata) delete mkpl_DDAdata;
  mkpl_DDAdata = new expdata_dda();


  mkpl_fm->GetExpDEtree()->SetBranchAddress("DE", &mkpl_DEdata);
  mkpl_fm->GetExpDAtree()->SetBranchAddress("DA", &mkpl_DAdata);
  mkpl_fm->GetExpDDtree()->SetBranchAddress("DD", &mkpl_DDdata);
  mkpl_fm->GetExpDDAtree()->SetBranchAddress("DDA", &mkpl_DDAdata);

  return;
}


TMultiGraph * mkplexpdatamanager::GetDE(const int theA, const int theZ, 
					const int ProjA, const int ProjZ, 
					const int TargA, const int TargZ, 
					const double ProjE, const bool LAB)
{
  TGraphErrors * exp_de = 0;
  TMultiGraph * result = new TMultiGraph();
  bool notitle(true);
  for (int i = 0; i < mkpl_fm->GetExpDEtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDEtree()->GetEntry(i);
      if (mkpl_DEdata->GetData()->GetProjectileZ() == ProjZ && 
	  mkpl_DEdata->GetData()->GetProjectileA() == ProjA && 
	  mkpl_DEdata->GetData()->GetTargetZ() == TargZ && 
	  mkpl_DEdata->GetData()->GetTargetA() == TargA && 
	  TMath::Abs(mkpl_DEdata->GetData()->GetProjectileEnergy()-ProjE) < 0.1)
	{
	  if (mkpl_DEdata->GetData()->GetParticleA() == theA &&
	      mkpl_DEdata->GetData()->GetParticleZ() == theZ && 
	      ((LAB || mkpl_DEdata->GetData()->IsCM()) &&  // XOR (if LAB=true that expr. will be
	       !(LAB && mkpl_DEdata->GetData()->IsCM()))   // true only if IsCM() return false)
	      )                                            
	    {
	      exp_de = (TGraphErrors*)(mkpl_DEdata->GetGraph());
	      result->Add(exp_de);
	      if (notitle)
		{
		  result->SetTitle(mkpl_DEdata->GetData()->GetGraphTitle());
		  notitle = false;
		}
	    }
	}
    }
  if (!result->GetListOfGraphs())
    {
      delete result;
      result = 0;
    }
  return result;
}




TMultiGraph * mkplexpdatamanager::GetDA(const int theA, const int theZ,
					const int ProjA, const int ProjZ, 
					const int TargA, const int TargZ, 
					const double ProjE, const bool LAB)
{
  TGraphErrors * exp_da = 0;
  TMultiGraph * result = new TMultiGraph();
  bool notitle(true);
  for (int i = 0; i < mkpl_fm->GetExpDAtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDAtree()->GetEntry(i);
      if (mkpl_DAdata->GetData()->GetProjectileZ() == ProjZ &&
	  mkpl_DAdata->GetData()->GetProjectileA() == ProjA &&
	  mkpl_DAdata->GetData()->GetTargetZ() == TargZ &&
	  mkpl_DAdata->GetData()->GetTargetA() == TargA &&
	  TMath::Abs(mkpl_DAdata->GetData()->GetProjectileEnergy()-ProjE)<0.1)
	{
	  if (mkpl_DAdata->GetData()->GetParticleA() == theA &&
	      mkpl_DAdata->GetData()->GetParticleZ() == theZ && 
	      ((LAB || mkpl_DAdata->GetData()->IsCM()) &&  // XOR 
	       !(LAB && mkpl_DAdata->GetData()->IsCM()))
	      )
	    {
	      exp_da = (TGraphErrors*)(mkpl_DAdata->GetGraph());
	      result->Add(exp_da);
	      if (notitle)
		{
		  result->SetTitle(mkpl_DAdata->GetData()->GetGraphTitle());
		  notitle = false;
		}
	    }
	}
    }
  if (!result->GetListOfGraphs())
    {
      delete result;
      result = 0;
    }
  return result;
}



vector<pair<double,TMultiGraph*> > * 
mkplexpdatamanager::GetDD(const int theA, const int theZ, 
			  const int ProjA, const int ProjZ, 
			  const int TargA, const int TargZ, 
			  const double ProjE, const bool LAB)
{
  TGraphErrors * exp_dd(0);
  map<double,TMultiGraph*> angles;
  for (int i = 0; i < mkpl_fm->GetExpDDtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDDtree()->GetEntry(i);
      if (mkpl_DDdata->GetData(0)->GetProjectileZ() == ProjZ &&
	  mkpl_DDdata->GetData(0)->GetProjectileA() == ProjA &&
	  mkpl_DDdata->GetData(0)->GetTargetZ() == TargZ &&
	  mkpl_DDdata->GetData(0)->GetTargetA() == TargA &&
	  TMath::Abs(mkpl_DDdata->GetData(0)->GetProjectileEnergy()-ProjE)<0.1)
	{
	  mkpl_fm->GetExpDDtree()->GetEntry(i);
	  if (mkpl_DDdata->GetData(0)->GetParticleA() == theA && 
	      mkpl_DDdata->GetData(0)->GetParticleZ() == theZ &&
	      ((LAB || mkpl_DDdata->GetData(0)->IsCM()) &&  // XOR
	       !(LAB && mkpl_DDdata->GetData(0)->IsCM()))
	      )
	    {
	      for (int a = 0; a < mkpl_DDdata->GetNangles(); a++)
		{
		  double ang = mkpl_DDdata->GetData(a)->GetAngle()*(TMath::Pi()/180.0);
		  map<double,TMultiGraph*>::iterator idx = angles.find(ang);
		  if (idx != angles.end())
		    {
		      exp_dd = (TGraphErrors*)(mkpl_DDdata->GetData(a)->GetGraph());
		      idx->second->Add(exp_dd);
		      exp_dd = 0;
		    }
		  else
		    {
		      exp_dd = (TGraphErrors*)(mkpl_DDdata->GetData(a)->GetGraph());
		      TMultiGraph * mg = new TMultiGraph();
		      mg->SetTitle(mkpl_DDdata->GetData(a)->GetGraphTitle());
		      mg->Add(exp_dd);
		      angles.insert(make_pair(ang,mg));
		      exp_dd = 0;
		    }
		}
	    }
	}
    }
  vector<pair<double,TMultiGraph*> > * result(0);
  if (!angles.empty())
    {
      result = new vector<pair<double,TMultiGraph*> >();
      for (map<double,TMultiGraph*>::iterator i = angles.begin();
	   i != angles.end(); ++i)
	{
	  result->push_back(make_pair(i->first,i->second));
	}
    }
  return result;
}


//vector<pair<pair<double,double>,TMultiGraph*> > * 
vector<mkplexpdda> * 
mkplexpdatamanager::GetDDA(const int theA, const int theZ, 
			   const int ProjA, const int ProjZ, 
			   const int TargA, const int TargZ, 
			   const double ProjE, const bool LAB)
{
  TGraphErrors * exp_dda(0);
  map<pair<double,double>,TMultiGraph*> ranges;
  for (int i = 0; i < mkpl_fm->GetExpDDAtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDDAtree()->GetEntry(i);
      if (mkpl_DDAdata->GetData(0)->GetProjectileZ() == ProjZ &&
	  mkpl_DDAdata->GetData(0)->GetProjectileA() == ProjA &&
	  mkpl_DDAdata->GetData(0)->GetTargetZ() == TargZ &&
	  mkpl_DDAdata->GetData(0)->GetTargetA() == TargA &&
	  TMath::Abs(mkpl_DDAdata->GetData(0)->GetProjectileEnergy()-ProjE)<0.1)
	{
	  mkpl_fm->GetExpDDAtree()->GetEntry(i);
	  if (mkpl_DDAdata->GetData(0)->GetParticleA() == theA && 
	      mkpl_DDAdata->GetData(0)->GetParticleZ() == theZ &&
	      ((LAB || mkpl_DDAdata->GetData(0)->IsCM()) &&  // XOR
	       !(LAB && mkpl_DDAdata->GetData(0)->IsCM()))
	      )
	    {
	      for (int r = 0; r < mkpl_DDAdata->GetNranges(); r++)
		{
		  double l = mkpl_DDAdata->GetData(r)->GetLowElimit();
		  double h = mkpl_DDAdata->GetData(r)->GetHighElimit();
		  pair<double,double> range = make_pair(l,h);
		  map<pair<double,double>,TMultiGraph*>::iterator idx = ranges.find(range);
		  if (idx != ranges.end())
		    {
		      exp_dda = (TGraphErrors*)(mkpl_DDAdata->GetData(r)->GetGraph());
		      idx->second->Add(exp_dda);
		      exp_dda = 0;
		    }
		  else
		    {
		      exp_dda = (TGraphErrors*)(mkpl_DDAdata->GetData(r)->GetGraph());
		      TMultiGraph * mg = new TMultiGraph();
		      mg->SetTitle(mkpl_DDAdata->GetData(r)->GetGraphTitle());
		      mg->Add(exp_dda);
		      ranges.insert(make_pair(range,mg));
		      exp_dda = 0;
		    }
		}
	    }
	}
    }
  //  vector<pair<pair<double,double>,TMultiGraph*> > * result(0);
  vector<mkplexpdda> * result(0);
  if (!ranges.empty())
    {
      //      result = new vector<pair<pair<double,double>,TMultiGraph*> >();
      result = new vector<mkplexpdda>();
      for (map<pair<double,double>,TMultiGraph*>::iterator i = ranges.begin();
	   i != ranges.end(); ++i)
	{
	  //	  result->push_back(make_pair(i->first,i->second));
	  mkplexpdda data(i->first.first, i->first.second, i->second);
	  result->push_back(data);
	}
    }
  return result;
}



TMultiGraph * mkplexpdatamanager::GetDD(const int energy)
{
  TMultiGraph * mg(0);
  if (energy < mkpl_fm->GetExpDDtree()->GetEntries())
    {
      mg = new TMultiGraph();
      mkpl_fm->GetExpDDtree()->GetEntry(energy);
      TString title(mkpl_DDdata->GetData(0)->GetGraphTitle());
      title.Remove(title.Index("MeV")+3);
      title += " #theta = ";
      char ang[5];
      int color(2);
      for (int angle = 0; angle < mkpl_DDdata->GetNangles(); angle++)
	{
	  TGraphErrors * graph = (TGraphErrors*)(mkpl_DDdata->GetData(angle)->GetGraph());
	  sprintf(ang,"%.1f",mkpl_DDdata->GetData(angle)->GetAngle());
	  title += TString(ang) + ", ";
	  graph->SetMarkerColor(color++);
	  graph->SetMarkerStyle(20);
	  graph->SetMarkerSize(0.75);
	  mg->Add(graph);
	}
      title.Remove(title.Last(','));
      mg->SetTitle(title);
    }
  return mg;
}


map<string,pair<int,int> > mkplexpdatamanager::WhichEjectiles(const double E)
{
  map<string,pair<int,int> >  theEjectiles;

  for (int i = 0; i < mkpl_fm->GetExpDEtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDEtree()->GetEntry(i);

      if (E <= 0.0 || TMath::Abs(E - mkpl_DEdata->GetData()->GetProjectileEnergy()) < 0.5)
	{
	  int A = mkpl_DEdata->GetData()->GetParticleA();
	  int Z = mkpl_DEdata->GetData()->GetParticleZ();
	  string symbol = mkpl_DEdata->GetData()->GetParticleSymbol().Data();

	  theEjectiles[symbol] = make_pair(A,Z);
	}
    }
  for (int i = 0; i < mkpl_fm->GetExpDAtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDAtree()->GetEntry(i);

      if (E <= 0.0 || TMath::Abs(E - mkpl_DAdata->GetData()->GetProjectileEnergy()) < 0.5)
	{
	  int A = mkpl_DAdata->GetData()->GetParticleA();
	  int Z = mkpl_DAdata->GetData()->GetParticleZ();
	  string symbol = mkpl_DAdata->GetData()->GetParticleSymbol().Data();
	  
	  theEjectiles[symbol] = make_pair(A,Z);
	}
    }
  for (int i = 0; i < mkpl_fm->GetExpDDtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDDtree()->GetEntry(i);
      
      if (E <= 0.0 || TMath::Abs(E - mkpl_DDdata->GetData(0)->GetProjectileEnergy()) < 0.5)
	{
	  int A = mkpl_DDdata->GetData(0)->GetParticleA();
	  int Z = mkpl_DDdata->GetData(0)->GetParticleZ();
	  string symbol = mkpl_DDdata->GetData(0)->GetParticleSymbol().Data();
	  
	  theEjectiles[symbol] = make_pair(A,Z);
	}
    }
  for (int i = 0; i < mkpl_fm->GetExpDDAtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDDAtree()->GetEntry(i);

      if (E <= 0.0 || TMath::Abs(E - mkpl_DDAdata->GetData(0)->GetProjectileEnergy()) < 0.5)
	{
	  int A = mkpl_DDAdata->GetData(0)->GetParticleA();
	  int Z = mkpl_DDAdata->GetData(0)->GetParticleZ();
	  string symbol = mkpl_DDAdata->GetData(0)->GetParticleSymbol().Data();
	  
	  theEjectiles[symbol] = make_pair(A,Z);
	}
    }

  return theEjectiles;

}


int mkplexpdatamanager::GetTargetZ() const
{
  int theZ(0);
  if (mkpl_fm->GetExpDEtree()->GetEntries() > 0)
    {
      mkpl_fm->GetExpDEtree()->GetEntry(0);
      theZ = mkpl_DEdata->GetData()->GetTargetZ();
    }
  else if (mkpl_fm->GetExpDAtree()->GetEntries() > 0)
    {
      mkpl_fm->GetExpDAtree()->GetEntry(0);
      theZ = mkpl_DAdata->GetData()->GetTargetZ();
    }
  else if (mkpl_fm->GetExpDDtree()->GetEntries() > 0)
    {
      mkpl_fm->GetExpDDtree()->GetEntry(0);
      theZ = mkpl_DDdata->GetHeader()->GetTargetZ();
    }
  else if (mkpl_fm->GetExpDDAtree()->GetEntries() > 0)
    {
      mkpl_fm->GetExpDDAtree()->GetEntry(0);
      theZ = mkpl_DDAdata->GetHeader()->GetTargetZ();
    }
  return theZ;
}


void mkplexpdatamanager::GetSummary(std::ostringstream & os)
{
  int TZ = this->GetTargetZ();
  if (mkpl_summary.empty()) this->Summarize();

  summary::iterator it;
  for (it = mkpl_summary.begin(); it != mkpl_summary.end(); ++it)
    {
      os << "Z = " << TZ << " A = " << it->first << "\n\tEnergies (MeV): ";
      set<double>::iterator last(it->second->end()); 
      --last;
      for (set<double>::iterator k = it->second->begin(); k != last; ++k)
	{
	  os << *k << ", ";
	}
      os <<  *last << '\n';
    }
  //  cout << os.rdbuf();
  return; 
}


void mkplexpdatamanager::Summarize()
{
  int nentries = int(mkpl_fm->GetExpDEtree()->GetEntries());
  for (int i = 0; i < nentries; i++)
    {
      mkpl_fm->GetExpDEtree()->GetEntry(i);
      int TA = mkpl_DEdata->GetHeader()->GetTargetA();
      double E = mkpl_DEdata->GetHeader()->GetProjectileEnergy();
      if (mkpl_summary.find(TA) == mkpl_summary.end())
	{
	  mkpl_summary.insert(make_pair(TA,new set<double>));
	}
      mkpl_summary[TA]->insert(E);
    }

  nentries = int(mkpl_fm->GetExpDAtree()->GetEntries());
  for (int i = 0; i < nentries; i++)
    {
      mkpl_fm->GetExpDAtree()->GetEntry(i);
      int TA = mkpl_DAdata->GetHeader()->GetTargetA();
      double E = mkpl_DAdata->GetHeader()->GetProjectileEnergy();
      if (mkpl_summary.find(TA) == mkpl_summary.end())
	{
	  mkpl_summary.insert(make_pair(TA,new set<double>));
	}
      mkpl_summary[TA]->insert(E);
    }


  nentries = int(mkpl_fm->GetExpDDtree()->GetEntries());
  for (int i = 0; i < nentries; i++)
    {
      mkpl_fm->GetExpDDtree()->GetEntry(i);
      int TA = mkpl_DDdata->GetHeader()->GetTargetA();
      double E = mkpl_DDdata->GetHeader()->GetProjectileEnergy();
      if (mkpl_summary.find(TA) == mkpl_summary.end())
	{
	  mkpl_summary.insert(make_pair(TA,new set<double>));
	}
      mkpl_summary[TA]->insert(E);
    }

  nentries = int(mkpl_fm->GetExpDDAtree()->GetEntries());
  for (int i = 0; i < nentries; i++)
    {
      mkpl_fm->GetExpDDAtree()->GetEntry(i);
      int TA = mkpl_DDAdata->GetHeader()->GetTargetA();
      double E = mkpl_DDAdata->GetHeader()->GetProjectileEnergy();
      if (mkpl_summary.find(TA) == mkpl_summary.end())
	{
	  mkpl_summary.insert(make_pair(TA,new set<double>));
	}
      mkpl_summary[TA]->insert(E);
    }
  
  return;
}


vector<double> mkplexpdatamanager::ListDEenergies()
{
  vector<double> result;
  for (int i = 0; i < mkpl_fm->GetExpDEtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDEtree()->GetEntry(i);
      result.push_back(mkpl_DEdata->GetHeader()->GetProjectileEnergy());
    }
  return result;
}

vector<double> mkplexpdatamanager::ListDAenergies()
{
  vector<double> result;
  for (int i = 0; i < mkpl_fm->GetExpDAtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDAtree()->GetEntry(i);
      result.push_back(mkpl_DAdata->GetHeader()->GetProjectileEnergy());
    }
  return result;
}

vector<double> mkplexpdatamanager::ListDDenergies()
{
  vector<double> result;
  for (int i = 0; i < mkpl_fm->GetExpDDtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDDtree()->GetEntry(i);
      result.push_back(mkpl_DDdata->GetHeader()->GetProjectileEnergy());
    }
  return result;
}

vector<double> mkplexpdatamanager::ListDDAenergies()
{
  vector<double> result;
  for (int i = 0; i < mkpl_fm->GetExpDDAtree()->GetEntries(); i++)
    {
      mkpl_fm->GetExpDDAtree()->GetEntry(i);
      result.push_back(mkpl_DDAdata->GetHeader()->GetProjectileEnergy());
    }
  return result;
}


vector<double> mkplexpdatamanager::ListDDangles(const int n)
{
  vector<double> result;
  if (n < mkpl_fm->GetExpDDtree()->GetEntries())
    {
      mkpl_fm->GetExpDDtree()->GetEntry(n);
    }
  for (int i = 0; i < mkpl_DDdata->GetNangles(); i++)
    {
      result.push_back(mkpl_DDdata->GetData(i)->GetAngle());
    }
  return result;
}


vector<pair<double,double> > mkplexpdatamanager::ListDDAranges(const int n)
{
  vector<pair<double,double> > result;
  if (n < mkpl_fm->GetExpDDAtree()->GetEntries())
    {
      mkpl_fm->GetExpDDAtree()->GetEntry(n);
    }
  for (int i = 0; i < mkpl_DDAdata->GetNranges(); i++)
    {
      result.push_back(make_pair<double,double>(mkpl_DDAdata->GetData(i)->GetLowElimit(),
						mkpl_DDAdata->GetData(i)->GetHighElimit()));
    }
  return result;
}



TGraphErrors * mkplexpdatamanager::GetDE(const int n)
{
  TGraphErrors * exp_de = 0;
  if (n < mkpl_fm->GetExpDEtree()->GetEntries())
    {
      mkpl_fm->GetExpDEtree()->GetEntry(n);
      exp_de = (TGraphErrors*)(mkpl_DEdata->GetGraph());
      exp_de->SetTitle(mkpl_DEdata->GetData()->GetGraphTitle());
    }
  return exp_de;
}


TGraphErrors * mkplexpdatamanager::GetDA(const int n)
{
  TGraphErrors * exp_da = 0;
  if (n < mkpl_fm->GetExpDAtree()->GetEntries())
    {
      mkpl_fm->GetExpDAtree()->GetEntry(n);
      exp_da = (TGraphErrors*)(mkpl_DAdata->GetGraph());
      exp_da->SetTitle(mkpl_DAdata->GetData()->GetGraphTitle());
    }
  return exp_da;
}

TGraphErrors * mkplexpdatamanager::GetDD(const int e, const int a)
{
  TGraphErrors * exp_dd = 0;
  if (e < mkpl_fm->GetExpDDtree()->GetEntries())
    {
      mkpl_fm->GetExpDDtree()->GetEntry(e);
      if (a < mkpl_DDdata->GetNangles())
	{
	  exp_dd = (TGraphErrors*)(mkpl_DDdata->GetData(a)->GetGraph());
	  TString title = mkpl_DDdata->GetData(a)->GetGraphTitle() +
	    "(" + mkpl_DDdata->GetData(a)->GetExforEntryCode() + ")";
	  exp_dd->SetTitle(title);
	}
    }
  return exp_dd;
}

TGraphErrors * mkplexpdatamanager::GetDDA(const int e, const int r)
{
  TGraphErrors * exp_dda = 0;
  if (e < mkpl_fm->GetExpDDAtree()->GetEntries())
    {
      mkpl_fm->GetExpDDAtree()->GetEntry(e);
      if (r < mkpl_DDAdata->GetNranges())
	{
	  exp_dda = (TGraphErrors*)(mkpl_DDAdata->GetData(r)->GetGraph());
	  TString title = mkpl_DDAdata->GetData(r)->GetGraphTitle() +
	    "(" + mkpl_DDAdata->GetData(r)->GetExforEntryCode() + ")";
	  exp_dda->SetTitle(title);
	}
    }
  return exp_dda;
}


void mkplexpdatamanager::GetDEdetails(const int n, std::ostringstream & os)
{
  if (n < mkpl_fm->GetExpDEtree()->GetEntries())
    {
      mkpl_fm->GetExpDEtree()->GetEntry(n);
      mkpl_DEdata->ShowYourSelf(os);
    }
  return;
}

void mkplexpdatamanager::GetDAdetails(const int n, std::ostringstream & os)
{
  if (n < mkpl_fm->GetExpDAtree()->GetEntries())
    {
      mkpl_fm->GetExpDAtree()->GetEntry(n);
      mkpl_DAdata->ShowYourSelf(os);
    }
  return;
}

void mkplexpdatamanager::GetDDdetails(const int e, const int a, std::ostringstream & os)
{
  if (e < mkpl_fm->GetExpDDtree()->GetEntries())
    {
      mkpl_fm->GetExpDDtree()->GetEntry(e);
      if (a < mkpl_DDdata->GetNangles())
	{
	  mkpl_DDdata->ShowYourSelf( a,os);
	}
    }
  return;
}



void mkplexpdatamanager::GetDDdetails(const int e, std::ostringstream & os)
{
  if (e < mkpl_fm->GetExpDDtree()->GetEntries())
    {
      mkpl_fm->GetExpDDtree()->GetEntry(e);
      mkpl_DDdata->ShowYourSelf(-1,os);
      for (int a = 0; a < mkpl_DDdata->GetNangles(); a++)
	{
	  os << "-----------------------------------------------\n";
	  mkpl_DDdata->ShowYourSelf( a,os);
	}
    }
  os << "-----------------------------------------------\n";
  return;
}


void mkplexpdatamanager::GetDDAdetails(const int e, const int r, std::ostringstream & os)
{
  if (e < mkpl_fm->GetExpDDAtree()->GetEntries())
    {
      mkpl_fm->GetExpDDAtree()->GetEntry(e);
      if (r < mkpl_DDAdata->GetNranges())
	{
	  mkpl_DDAdata->ShowYourSelf( r,os);
	}
    }
  return;
}

void mkplexpdatamanager::GetDDAdetails(const int e, std::ostringstream & os)
{
  if (e < mkpl_fm->GetExpDDAtree()->GetEntries())
    {
      mkpl_fm->GetExpDDAtree()->GetEntry(e);
      mkpl_DDAdata->ShowYourSelf(-1,os);
      for (int r = 0; r < mkpl_DDAdata->GetNranges(); r++)
	{
	  os << "-----------------------------------------------\n";
	  mkpl_DDAdata->ShowYourSelf( r,os);
	}
    }
  os << "-----------------------------------------------\n";
  return;
}
