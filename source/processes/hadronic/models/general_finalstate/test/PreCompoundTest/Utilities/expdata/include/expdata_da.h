#ifndef expdata_da_h
#define expdata_da_h

#include "expdata.h"
#include "expdata_daunit.h"
#include "TObject.h"
#include "TGraph.h"
#include <iostream>
#include <sstream>

class expdata_da : public TObject
{ 
public:
  // Default constructor
  inline expdata_da();
  
  // Destructor
  inline ~expdata_da(); 
  
  inline void Reset();

  inline expdata_daunit * GetData() const;

  inline void SetData(const Int_t n,
		      const Double_t * e, const Double_t * ee,
		      const Double_t * d, const Double_t * ed);
    
  inline void SetData(const Int_t n,
		      const Double_t * e, const Double_t * eel, const Double_t * eeh,
		      const Double_t * d, const Double_t * edl, const Double_t * edh);

  const TGraph * GetGraph() const;
    
  inline expdata * GetHeader();

  inline void ShowYourSelf() const;
  void ShowYourSelf(std::ostringstream & os) const;

private:

  expdata     theDataHeader;
  expdata_daunit *      dat;               //->

  ClassDef(expdata_da,1)
};

#include "expdata_da.icc"

#endif

