#ifndef mkplsimfilemanager_h
#define mkplsimfilemanager_h

#include <string>

using namespace std;

#include "TFile.h"
#include "TTree.h"

class mkplsimfilemanager
{
public:
  mkplsimfilemanager(const string & simfilename, const int level=0);

  inline ~mkplsimfilemanager();

  inline TTree * GetTree();

  inline string GetFileName();
  inline string GetBaseName();

  inline void SetVerbosity(const int level);
  inline int  GetVerbosity() const;

private:
  inline mkplsimfilemanager();

  inline mkplsimfilemanager(const mkplsimfilemanager& right);

  inline const mkplsimfilemanager& operator=(const mkplsimfilemanager& right);

private:
  TFile * mkpl_simfile;
  
  TTree * mkpl_simtree;

  string  mkpl_basename;
  int     mkpl_verbose;
};

#include "mkplsimfilemanager.icc"

#endif
