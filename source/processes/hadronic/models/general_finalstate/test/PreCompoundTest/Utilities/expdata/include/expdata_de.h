#ifndef expdata_de_h
#define expdata_de_h

#include "expdata.h"
#include "expdata_deunit.h"
#include "TObject.h"
#include "TGraph.h"
#include <iostream>
#include <sstream>

class expdata_de : public TObject
{ 
public:
  // Default constructor
  inline expdata_de();
  
  // Destructor
  inline ~expdata_de();

  inline void Reset();

  inline expdata_deunit * GetData() const;

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
  
  expdata               theDataHeader;
  expdata_deunit *      dat;               //->

  ClassDef(expdata_de,1)
};

#include "expdata_de.icc"

#endif

