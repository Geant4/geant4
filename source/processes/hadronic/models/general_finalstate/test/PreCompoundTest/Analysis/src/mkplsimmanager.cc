#include "mkplsimmanager.h"
#include "TranslateFragment.h"

#include <iostream>
#include <map>

#include "TVector3.h"
#include "TMultiGraph.h"
#include "TPostScript.h"
#include "TApplication.h"



void mkplsimmanager::DeleteHistograms()
{
  if (!mkpl_comparisons.empty()) 
    for_each(mkpl_comparisons.begin(),mkpl_comparisons.end(),DeleteComparison());
  if (mkpl_tests) delete mkpl_tests;
  return;
}



void mkplsimmanager::Initialize(mkplfilemanager * fm)
{
  // Delete Histograms if they already exists
  this->DeleteHistograms();

  // Get the simulation tree
  mkpl_fm = fm;
  if (mkpl_reaction) delete mkpl_reaction;
  mkpl_reaction = new precoreaction();
  mkpl_fm->GetSimtree()->SetBranchAddress("reaction",&mkpl_reaction);
    
  mkpl_fm->GetSimtree()->GetEntry(0);
  // Get info for experimental comparison
  if (mkpl_reaction->IsNatural())
    mkpl_target_A = 0;
  else 
    mkpl_target_A = mkpl_reaction->GetTargetA();
  mkpl_target_Z = mkpl_reaction->GetTargetZ();
  mkpl_projectile_A = mkpl_reaction->GetProjectileA();
  mkpl_projectile_Z = mkpl_reaction->GetProjectileZ(); 
  mkpl_E = mkpl_reaction->GetProjectileKineticE();
  return;
}


void mkplsimmanager::PrepareTestHistograms()
{
  if (mkpl_tests) delete mkpl_tests;
  mkpl_fm->GetSimtree()->GetEntry(0);
  TVector3 IncidentDirection = mkpl_reaction->GetProjectileMomentum().Vect().Unit();
  // run a first loop in order to get the particle and process types
  set<string> particletypes;
  set<string> processtypes;
  for (int i = 0; i < mkpl_fm->GetSimtree()->GetEntries(); i++)
    {
      mkpl_fm->GetSimtree()->GetEntry(i);
      for (int j = 0; j < mkpl_reaction->GetNumOfFragments(); j++)
	{
	  precofragment * fragment = mkpl_reaction->GetFragment(j);
	  particletypes.insert(TranslateFragment(fragment->GetFragmentName()));
	  processtypes.insert(fragment->GetProcessName());
	}
    }

  mkpl_tests = new mkpltesthistograms();
  mkpl_tests->PrepareHistograms(IncidentDirection,particletypes,processtypes);
  return;
}
    
void mkplsimmanager::PrepareComparisonHistograms(mkplexpdatamanager * expdata, bool components)
{
  // Get info for histograms preparation
  mkpl_fm->GetSimtree()->GetEntry(0);
  // Cross Section in millibarn
  double xs = 1.e3*mkpl_reaction->GetCrossSection();
  // Number of simulated events
  double entries = mkpl_fm->GetSimtree()->GetEntries();

  // Find out which kind of ejectiles there are in experimental file
  map<string,pair<int,int> > ejectiles = expdata->WhichEjectiles(mkpl_E);
  // Delete (if they exists) comparison histograms
  if (!mkpl_comparisons.empty()) 
    for_each(mkpl_comparisons.begin(),mkpl_comparisons.end(),DeleteComparison());

  // Create new Histograms for Comparisons
  map<string,pair<int,int> >::iterator i;
  for (i=ejectiles.begin(); i!=ejectiles.end();++i)
    {
      mkplcomparisonhistograms * comparison = new mkplcomparisonhistograms(i->second.first,i->second.second,i->first.c_str());
      comparison->SetProjectileA(mkpl_projectile_A);
      comparison->SetProjectileZ(mkpl_projectile_Z);
      comparison->SetTargetA(mkpl_target_A);
      comparison->SetTargetZ(mkpl_target_Z);
      comparison->SetReactionE(mkpl_E);
      comparison->PrepareHistograms(expdata,components,xs,entries);
      mkpl_comparisons.push_back(comparison);
    }

  return;
}


void mkplsimmanager::FillHistograms()
{
  // Loop for simulated events
  for (int i = 0; i < mkpl_fm->GetSimtree()->GetEntries(); i++)
    {
      mkpl_fm->GetSimtree()->GetEntry(i);

      // Loop for all fragments in the i-th reaction
      for (int j = 0; j < mkpl_reaction->GetNumOfFragments(); j++)
	{
	  precofragment * fragment = mkpl_reaction->GetFragment(j);

	  TVector3 IncidentDirection = mkpl_reaction->GetProjectileMomentum().Vect().Unit();
	  TVector3 boostToCM = mkpl_reaction->GetBoostToCM();
	  double theE = fragment->GetKineticE();
	  double theEinCM = fragment->GetKineticE(boostToCM);
	  double Theta = fragment->GetTheta(IncidentDirection);
	  double ThetaCM = fragment->GetTheta(IncidentDirection,boostToCM);
	  int theA = fragment->GetA();
	  int theZ = fragment->GetZ();
	  string procname = fragment->GetProcessName();

	  // Fill the comparison histograms
	  // ------------------------------
	  for (vector<mkplcomparisonhistograms*>::iterator i = mkpl_comparisons.begin();
	       i != mkpl_comparisons.end(); ++i)
	    (*i)->Fill(theA,theZ,theE,Theta,theEinCM,ThetaCM,procname);

	  if (mkpl_tests)
	    {
	      string fname = TranslateFragment(fragment->GetFragmentName());
	      TVector3 fragP = fragment->GetMomentum().Vect();
	      
	      // Fill the test histograms
	      //-------------------------
	      mkpl_tests->FillHistograms(Theta,fragP,procname,fname);
	    }
	}
      if (mkpl_tests)
	{
	  int deltaA = mkpl_reaction->CheckConservationA();
	  int deltaZ = mkpl_reaction->CheckConservationZ();
	  double deltaE = mkpl_reaction->CheckConservationE();
	  double deltaP = mkpl_reaction->CheckConservationPMag();
	  TVector3 MomentumTest = mkpl_reaction->CheckConservationP();
	  // Fill the conservation test histograms
	  // -------------------------------------
	  mkpl_tests->FillConservationHistograms(deltaA,deltaZ,deltaE,deltaP,MomentumTest);
	}
    }
  return;
}
