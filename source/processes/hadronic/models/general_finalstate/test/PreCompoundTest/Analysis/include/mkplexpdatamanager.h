#ifndef mkplexpdatamanager_h
#define mkplexpdatamanager_h

#include "mkplfilemanager.h"
#include "mkplexpdda.h"
#include "expdata_de.h"
#include "expdata_da.h"
#include "expdata_dd.h"
#include "expdata_dda.h"

#include "TMultiGraph.h"

#include <utility>
#include <vector>
#include <map>
#include <set>
#include <sstream>

using namespace std;

class mkplexpdatamanager
{
public:
  inline mkplexpdatamanager();
  inline mkplexpdatamanager(mkplfilemanager & fm);
  inline mkplexpdatamanager(mkplfilemanager * fm);
  inline ~mkplexpdatamanager();

  void Initialize(mkplfilemanager * fm);

  TMultiGraph * GetDE(const int theA, const int theZ, 
		      const int ProjA, const int ProjZ, 
		      const int TargA, const int TargZ, 
		      const double ProjE, const bool LAB=true);
  
  TMultiGraph * GetDA(const int theA, const int theZ, 
		      const int ProjA, const int ProjZ, 
		      const int TargA, const int TargZ, 
		      const double ProjE, const bool LAB=true);
  
  vector<pair<double,TMultiGraph*> > * GetDD(const int theA, const int theZ, 
					     const int ProjA, const int ProjZ, 
					     const int TargA, const int TargZ, 
					     const double ProjE, const bool LAB=true);

  vector<mkplexpdda> * GetDDA(const int theA, const int theZ, 
			      const int ProjA, const int ProjZ, 
			      const int TargA, const int TargZ, 
			      const double ProjE, const bool LAB=true);

  //  vector<pair<pair<double,double>,TMultiGraph*> > * GetDDA(const int theA, const int theZ, 
  //						   const int ProjA, const int ProjZ, 
  //						   const int TargA, const int TargZ, 
  //						   const double ProjE, const bool LAB=true);

  


  TGraphErrors * GetDE(const int n);
  TGraphErrors * GetDA(const int n);
  TGraphErrors * GetDD(const int e, const int a);
  TMultiGraph  * GetDD(const int energy);
  TGraphErrors * GetDDA(const int e, const int r);

  int GetTargetZ() const;
  map<string,pair<int,int> > WhichEjectiles(const double E=-1.0);
  void GetSummary(std::ostringstream & os);
  vector<double> ListDEenergies();
  vector<double> ListDAenergies();
  vector<double> ListDDenergies();
  vector<double> ListDDAenergies();
  vector<double> ListDDangles(const int n);
  vector<pair<double,double> > ListDDAranges(const int n);
  inline mkplfilemanager * GetFileManager() const;
  void GetDEdetails(const int n, std::ostringstream& os);
  void GetDAdetails(const int n, std::ostringstream& os);
  void GetDDdetails(const int e, std::ostringstream& os);
  void GetDDdetails(const int e, const int a, std::ostringstream& os);
  void GetDDAdetails(const int e, std::ostringstream& os);
  void GetDDAdetails(const int e, const int a, std::ostringstream& os);

private:
  inline mkplexpdatamanager(const mkplexpdatamanager& right);

  inline const mkplexpdatamanager& operator=(const mkplexpdatamanager& right);

  void Summarize();

private:
  typedef map<int, set<double>* > summary;

  mkplfilemanager * mkpl_fm;
  expdata_de      * mkpl_DEdata;
  expdata_da      * mkpl_DAdata;
  expdata_dd      * mkpl_DDdata; 
  expdata_dda     * mkpl_DDAdata; 
  
  summary           mkpl_summary;

};

#include "mkplexpdatamanager.icc"

#endif
