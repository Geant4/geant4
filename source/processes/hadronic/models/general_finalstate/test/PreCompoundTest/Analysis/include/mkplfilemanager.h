#ifndef mkplfilemanager_h 
#define mkplfilemanager_h 

#include <string>

using namespace std;

#include "TFile.h"
#include "TTree.h"

#include "mkplexpfilemanager.h"
#include "mkplsimfilemanager.h"

class mkplfilemanager
{
public:
  inline mkplfilemanager();
  inline ~mkplfilemanager();

  inline void OpenSimFile(const string & filename);
  inline void OpenExpFile(const string & filename);

  inline void CloseSimFile();
  inline void CloseExpFile(); 

  inline TTree * GetExpDEtree();
  inline TTree * GetExpDAtree();
  inline TTree * GetExpDDtree();
  inline TTree * GetExpDDAtree();
  inline TTree * GetSimtree();

  inline string GetSimFileName();
  inline string GetExpFileName();

  inline mkplexpfilemanager * GetExpFile();
  inline mkplsimfilemanager * GetSimFile();

  inline void SetVerbosity(const int level);

private:
  inline mkplfilemanager(const mkplfilemanager& right);
  inline const mkplfilemanager& operator=(const mkplfilemanager& right);

private:
    
    mkplexpfilemanager * mkpl_expfile;
    mkplsimfilemanager * mkpl_simfile;

};

#include "mkplfilemanager.icc"

#endif
