#ifndef mkplexpfilemanager_h
#define mkplexpfilemanager_h

#include <string>

using namespace std;


#include "TFile.h"
#include "TTree.h"

class mkplexpfilemanager
{
public:

  mkplexpfilemanager(const string& expfile, const int level=0);
  inline ~mkplexpfilemanager();

  inline TTree * GetDEtree();
  inline TTree * GetDAtree();
  inline TTree * GetDDtree();
  inline TTree * GetDDAtree();

  inline string GetFileName();
  inline string GetBaseName();

  inline void SetVerbosity(const int level);
  inline int  GetVerbosity() const;

private:
  inline mkplexpfilemanager();

  inline mkplexpfilemanager(const mkplexpfilemanager& right);

  inline const mkplexpfilemanager& operator=(const mkplexpfilemanager& right); 

private:
  TFile * mkpl_expfile;

  TTree * mkpl_exptree_de;
  TTree * mkpl_exptree_da;
  TTree * mkpl_exptree_dd;
  TTree * mkpl_exptree_dda;

  string  mkpl_basename;
  int     mkpl_verbose;
};

#include "mkplexpfilemanager.icc"

#endif
