#ifndef expdata_daunit_h
#define expdata_daunit_h

#include "expdata_unit.h"
#include "TGraphErrors.h"
#include <iostream>
#include <sstream>

class expdata_daunit : public expdata_unit
{
public:
  // Default constructor
  inline expdata_daunit();
  // Copy constructor
  inline expdata_daunit(const expdata_daunit & val);
  // Destructor
  inline ~expdata_daunit();
  // Assignment Operator 
  inline const expdata_daunit & operator=(const expdata_daunit & val);
  
  inline expdata_daunit(const Int_t n, 
			const Double_t * e, 
			const Double_t * ee,
			const Double_t * d, 
			const Double_t * ed);

  
  inline expdata_daunit(const Int_t n, 
			const Double_t * e, 
			const Double_t * eel, 
			const Double_t * eeh,
			const Double_t * d, 
			const Double_t * edl, 
			const Double_t * edh);
  
  
  inline Int_t GetId() const;
  
  inline Double_t GetAngle(const Int_t i) const;
  inline Double_t GetEngleError(const Int_t i) const;
  
  inline Double_t GetCrossSection(const Int_t i) const;
  inline Double_t GetCrossSectionError(const Int_t i) const;
  
  TString GetGraphTitle() const;

  inline void ShowYourSelf() const;

  void ShowYourSelf(std::ostringstream & os) const;
  
private:
  Int_t Id;
  
  ClassDef(expdata_daunit,1)
};

#include "expdata_daunit.icc"

#endif



