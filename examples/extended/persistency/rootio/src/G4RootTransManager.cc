// $Id: G4RootTransManager.cc,v 1.1 2002-12-04 02:44:29 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4RootTransManager.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4RootTransManager.hh"

// Addtional Include:
#include "G4RootIOManager.hh"
#include "G4PersistencyCenter.hh"
#include "TFile.h"
#include "TTree.h"

// Implementation of Constructor #1
G4RootTransManager::G4RootTransManager()
 : f_lastreadfile(0), f_lastwritefile(0), f_lastreadtree(0), f_lastwritetree(0), f_treeName("FADS"), f_treeDesc("FADS output tree"), f_bufsize(64000), f_split(1), f_comp(1)
{
  Initialize();
}

// Implementation of SelectReadFile
bool G4RootTransManager::SelectReadFile(std::string obj, std::string file)
{
  // obtain TFile pointer for "file" for reading
  TFile* tfile = LookUp(file, "READ");
  if ( tfile == 0 ) return false;

  // cd to the file
  tfile->cd();

  // Get the tree "FADS"
  TTree* tree = (TTree*)tfile->Get(f_treeName.c_str());
  if ( tree == 0 ) return false;

  // cache the pointer to the last write file
  f_lastreadfile = tfile;
  f_lastreadtree = tree;

  return true;
}

// Implementation of SelectWriteFile
bool G4RootTransManager::SelectWriteFile(std::string obj, std::string file)
{
  // obtain TFile pointer for "obj" and "file" for writing
  TFile* tfile = LookUp(file, "WRITE");
  if ( tfile == 0 ) return false;

  // cd to the file
  tfile->cd();

  // Get the tree "FADS"
  TTree* tree = (TTree*)tfile->Get(f_treeName.c_str());
  if ( tree == 0 ) {
    tree = new TTree(f_treeName.c_str(),f_treeDesc.c_str());
    if ( tree == 0 ) return false;
    tree->SetAutoSave(1000000000);  // autosave when 1 Gbyte written
    if ( m_verbose > 2 )
      std::cout << "G4RootTransManager: TTree \"" << f_treeName << "\" created."
                << std::endl;
  }

  // cache the pointer to the last write file
  f_lastwritefile = tfile;
  f_lastwritetree = tree;

  return true;
}

// Implementation of Commit
void G4RootTransManager::Commit()
{
  if( f_lastwritefile != 0 ) {
    // f_lastwritefile->cd();
    // if( f_lastwritetree != 0 ) f_lastwritetree->Fill();
    f_lastwritefile->Write();
  }
}

// Implementation of Initialize
void G4RootTransManager::Initialize()
{
  // Initialize the root file catalog
  RootFileCatalog::iterator itr;
  for ( itr = f_catalog.begin(); itr != f_catalog.end(); itr++ ) {
    (*itr).second = 0;
  }
}

// Implementation of BufferSize
int G4RootTransManager::BufferSize()
{
  if (f_split) return f_bufsize/4;
  return f_bufsize;
}

// Implementation of LookUp
TFile* G4RootTransManager::LookUp(string file, string mode)
{
  if ( (*(f_catalog.find(file))).second != 0 ) {
    return f_catalog[file];
  } else {
    if ( mode == "WRITE" ) {
      f_catalog[file] = new TFile(file.c_str(), "UPDATE");

      cout << "LookUp: File: " << file << " for write" << endl;
      f_catalog[file]->SetCompressionLevel(f_comp);
      return f_catalog[file];
    } else if ( mode == "READ" ) {
      f_catalog[file] = new TFile(file.c_str(), "READ");

      cout << "LookUp: File: " << file << " for read" << endl;
      return f_catalog[file];
    } else {
      cout << "LookUp: File: " << file << ", f_catalog=0" << endl;
      return 0;
    }
  }
}

// End of G4RootTransManager.cc

