#include "mkplsimfilemanager.h"

#include <iostream>

mkplsimfilemanager::mkplsimfilemanager(const string & simfilename, const int level)
  : mkpl_simfile(0), mkpl_simtree(0), mkpl_verbose(level)
{
  //++++++++++++++++++++++++++++
  //+ Open the simulation tree +
  //++++++++++++++++++++++++++++
  // Open the simulation data file 
  mkpl_simfile = new TFile(simfilename.c_str(),"READ");
  if (!mkpl_simfile) 
    {
      cerr << "Can not open simulated data file: " << simfilename << '\n';
      exit(2);
    }

  // Get the simulated tree
  if (mkpl_verbose > 0) cout << "Loading simulated tree...";
  mkpl_simtree = (TTree*)(mkpl_simfile->Get("PCT"));
  if (mkpl_verbose > 0) cout << " done!\n"; 

  // The base name
  mkpl_basename = simfilename;
  string::size_type idx = mkpl_basename.find(".root");
  if (idx != string::npos) 
    {
      mkpl_basename.replace(idx,5,"");
    }
  idx = mkpl_basename.find_last_of('/');
  if (idx != string::npos) 
    {
      mkpl_basename.replace(0,idx+1,"");
    }

}
