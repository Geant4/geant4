#ifndef expdata_dd_h
#define expdata_dd_h

#include "expdata.h"
#include "expdata_ddunit.h"
#include "TObject.h"
#include "TGraph.h"
#include "TObjArray.h"
#include <iostream>
#include <sstream>

class expdata_dd : public TObject
{
public:
  // Default Constructor
  expdata_dd();

  // Destructor
  inline ~expdata_dd();
  
  void Reset();
  
  Int_t GetNangles() const;
  
  TObjArray * GetData() const;
  
  expdata_ddunit * GetData(const Int_t i) const;
  expdata_ddunit * GetData(const Double_t ang) const;
  
  inline void AddAngle(const Int_t n, const Double_t angle, 
		       const Double_t * e, const Double_t * ee,
		       const Double_t * d, const Double_t * ed);
  
  inline void AddAngle(const Int_t n, const Double_t angle, 
		       const Double_t * e, const Double_t * eel, const Double_t * eeh,
		       const Double_t * d, const Double_t * edl, const Double_t * edh);
  
  const TGraph * GetGraph(const Int_t i) const;
  
  inline expdata * GetHeader();

  inline void ShowYourSelf(const Int_t v) const;
  void ShowYourSelf(const Int_t v, std::ostringstream & os) const;

private:
  Int_t                 nangles;
  expdata               theDataHeader;
  TObjArray *           dat;               //->
  
  ClassDef(expdata_dd,1)
};

#include "expdata_dd.icc"

#endif





