#ifndef expdata_deunit_h
#define expdata_deunit_h

#include "expdata_unit.h"
#include "TGraphErrors.h"
#include <iostream>
#include <sstream>

class expdata_deunit : public expdata_unit
{
public:
  // Default constructor
  inline expdata_deunit();
  // Copy constructor
  inline expdata_deunit(const expdata_deunit & val);
  // Destructor
  inline ~expdata_deunit();
  // Assignment Operator 
  inline const expdata_deunit & operator=(const expdata_deunit & val);
  
  inline expdata_deunit(const Int_t n, 
			const Double_t * e, const Double_t * ee,
			const Double_t * d, const Double_t * ed);

  
  inline expdata_deunit(const Int_t n, 
			const Double_t * e, const Double_t * eel, const Double_t * eeh,
			const Double_t * d, const Double_t * edl, const Double_t * edh);
  
  inline Int_t GetId() const;
  
  inline Double_t GetEnergy(const Int_t i) const;
  inline Double_t GetEnergyError(const Int_t i) const;
  
  inline Double_t GetCrossSection(const Int_t i) const;
  inline Double_t GetCrossSectionError(const Int_t i) const;
  
  TString GetGraphTitle() const;

  inline void ShowYourSelf() const;
  void ShowYourSelf(std::ostringstream & os) const;
  
private:
  Int_t Id;
  
  ClassDef(expdata_deunit,1)
};

#include "expdata_deunit.icc"

#endif



