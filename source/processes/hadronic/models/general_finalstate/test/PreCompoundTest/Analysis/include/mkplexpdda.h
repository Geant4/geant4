#ifndef mkplexpdda_h
#define mkplexpdda_h

#include "TMultiGraph.h"

using namespace std;

class mkplexpdda
{
public:
  inline mkplexpdda();
  
  inline ~mkplexpdda();
  
  inline mkplexpdda(const double li, const double ls, TMultiGraph * ptr);

  inline mkplexpdda(const mkplexpdda& right);

  inline const mkplexpdda& operator=(const mkplexpdda & right);

  inline double GetInfCut() const;

  inline double GetSupCut() const;

  inline void SetRange(const double li, const double ls);

  inline void SetLimits(const double li, const double ls);

  inline TMultiGraph * GetData();

  inline void SetData(TMultiGraph * ptr);

private:

  double mkpl_liminf;
  double mkpl_limsup;
  TMultiGraph * mkpl_mg;
};

#include "mkplexpdda.icc"

#endif
