#include "mkplfilemanager.h"

#include <iostream>


mkplfilemanager::mkplfilemanager(const string& expfilename, const string& simfilename)
    : mkpl_expfile(0), mkpl_simfile(0), mkpl_exptree_de(0), mkpl_exptree_da(0),
      mkpl_exptree_dd(0), mkpl_simtree(0)
{
    //+++++++++++++++++++++++++++++++++++
    //+ Open the experimental data file +
    //+++++++++++++++++++++++++++++++++++
    // open the experimental file 
    mkpl_expfile = new TFile(expfilename.c_str());
    if (!mkpl_expfile)
    {
	cerr << "Can not open experimental data file: " << expfilename << '\n';
	exit(2);
    }

    // generate the tree names
    string basexptree = mkpl_expfile->GetName();
    string::size_type idx = basexptree.find(".root");
    if (idx != string::npos) 
    {
	basexptree.replace(idx,5,"");
    }
    idx = basexptree.find_last_of('/');
    if (idx != string::npos) 
    {
	basexptree.replace(0,idx+1,"");
    }
    
    string tDEname = basexptree + "_DE";
    string tDAname = basexptree + "_DA";
    string tDDname = basexptree + "_DD";

    // Experimental trees
    cout << "Loading experimental trees...";
    mkpl_exptree_de = (TTree*)(mkpl_expfile->Get(tDEname.c_str()));
    mkpl_exptree_da = (TTree*)(mkpl_expfile->Get(tDAname.c_str()));
    mkpl_exptree_dd = (TTree*)(mkpl_expfile->Get(tDDname.c_str()));
    cout << " done!\n";

    //++++++++++++++++++++++++++++
    //+ Open the simulation tree +
    //++++++++++++++++++++++++++++
    // Open the simulation data file 
    mkpl_simfile = new TFile(simfilename.c_str());
    if (!mkpl_simfile) 
    {
	cerr << "Can not open simulated data file: " << simfilename << '\n';
	exit(3);
    }

    // Get the simulated tree
    cout << "Loading simulated tree...";
    mkpl_simtree = (TTree*)(mkpl_simfile->Get("PCT"));
    cout << " done!\n";
}

mkplfilemanager::~mkplfilemanager()
{
    if (mkpl_expfile) delete mkpl_expfile;
    if (mkpl_simfile) delete mkpl_simfile;
}
