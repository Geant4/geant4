#ifndef expdata_dda_h
#define expdata_dda_h

#include "expdata.h"
#include "expdata_ddaunit.h"
#include "TObject.h"
#include "TGraph.h"
#include "TObjArray.h"
#include <iostream>
#include <sstream>

class expdata_dda : public TObject
{
public:
  // Default Constructor
  expdata_dda();
  
  // Destructor 
  inline ~expdata_dda();

  void Reset();

  Int_t GetNranges() const;
  
  TObjArray * GetData() const;

  expdata_ddaunit * GetData(const Int_t i) const;
  expdata_ddaunit * GetData(const Double_t emin, const Double_t emax);

  inline void AddRange(const Int_t n, const Double_t Emin, const Double_t Emax,
		       const Double_t * a, const Double_t * ea,
		       const Double_t * d, const Double_t * ed);

  inline void AddRange(const Int_t n, const Double_t Emin, const Double_t Emax,
		       const Double_t * a, const Double_t * eal, const Double_t * eah,
		       const Double_t * d, const Double_t * edl, const Double_t * edh);

  const TGraph * GetGraph(const Int_t i) const;

  inline expdata * GetHeader();

  inline void ShowYourSelf(const Int_t v) const;
  void ShowYourSelf(const Int_t v, std::ostringstream & os) const;
  
private:
  Int_t       nranges;
  expdata     theDataHeader;
  TObjArray * dat;                //->

  ClassDef(expdata_dda,1)
};

#include "expdata_dda.icc"

#endif
