#ifndef mkplsimmanager_h
#define mkplsimmanager_h

#include <vector>
#include <set>
#include <string>
#include <utility>
#include <sstream>

#include "mkplfilemanager.h"
#include "mkplexpdatamanager.h"
#include "mkpltesthistograms.h"
#include "mkplcomparisonhistograms.h"
#include "precoreaction.h"


class TApplication;

using namespace std;

class mkplsimmanager 
{
public:
  inline mkplsimmanager();    
  inline mkplsimmanager(mkplfilemanager & fm);
  inline mkplsimmanager(mkplfilemanager * fm);
  inline ~mkplsimmanager();
  
  void Initialize(mkplfilemanager * fm);
  
  void PrepareTestHistograms();

  void PrepareComparisonHistograms(mkplexpdatamanager * expdata, bool components = true);
  
  void FillHistograms();
  
  inline void PrintInfo() const; 
  inline void GetSummary(std::ostringstream & os) const;

  inline int GetTargetZ() const;
  inline mkpltesthistograms * GetTestHistograms();
  inline vector<mkplcomparisonhistograms*> * GetComparisonHistograms();


private:
  inline mkplsimmanager(const mkplsimmanager& right);
  
  inline const mkplsimmanager& operator=(const mkplsimmanager& right);
  

  void DeleteHistograms();
      
private:
  mkplfilemanager * mkpl_fm;
  precoreaction * mkpl_reaction;
  int mkpl_projectile_A;
  int mkpl_projectile_Z;
  int mkpl_target_A;
  int mkpl_target_Z;
  double mkpl_E;
  
  // Basic Test Histograms
  // ---------------------
  mkpltesthistograms * mkpl_tests;
  
  // Comparison with Experimental Data Histograms
  // --------------------------------------------
  vector<mkplcomparisonhistograms*> mkpl_comparisons;

  struct DeleteComparison
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };

};

#include "mkplsimmanager.icc"

#endif
