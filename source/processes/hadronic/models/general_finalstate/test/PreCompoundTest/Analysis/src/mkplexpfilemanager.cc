#include "mkplexpfilemanager.h"

#include <iostream>

mkplexpfilemanager::mkplexpfilemanager(const string& expfilename, const int level)
  : mkpl_expfile(0), mkpl_exptree_de(0), mkpl_exptree_da(0), mkpl_exptree_dd(0), 
    mkpl_exptree_dda(0), mkpl_verbose(level)
{
  //+++++++++++++++++++++++++++++++++++
  //+ Open the experimental data file +
  //+++++++++++++++++++++++++++++++++++
  // open the experimental file 
  mkpl_expfile = new TFile(expfilename.c_str()),"READ";
  if (!mkpl_expfile)
    {
      cerr << "Can not open experimental data file: " << expfilename << '\n';
      exit(2);
    }

  
  // generate the tree names
  mkpl_basename = mkpl_expfile->GetName();
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


  string tDEname = mkpl_basename + "_DE";
  string tDAname = mkpl_basename + "_DA";
  string tDDname = mkpl_basename + "_DD";
  string tDDAname = mkpl_basename + "_DDA";

  // Experimental trees
  if (mkpl_verbose > 0) cout << "Loading experimental trees...";
  mkpl_exptree_de = (TTree*)(mkpl_expfile->Get(tDEname.c_str()));
  mkpl_exptree_da = (TTree*)(mkpl_expfile->Get(tDAname.c_str()));
  mkpl_exptree_dd = (TTree*)(mkpl_expfile->Get(tDDname.c_str()));
  mkpl_exptree_dda = (TTree*)(mkpl_expfile->Get(tDDAname.c_str()));
  if (mkpl_verbose) cout << " done!\n";

}



