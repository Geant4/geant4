#ifndef expdata_ddunit_h
#define expdata_ddunit_h

#include "expdata_unit.h"
#include "TGraphErrors.h"
#include <iostream>
#include <sstream>

class expdata_ddunit : public expdata_unit
{
public:
  // Default constructor
  inline expdata_ddunit();
  // Copy constructor
  inline expdata_ddunit(const expdata_ddunit & val);
  // Destructor
  inline ~expdata_ddunit();
  // Assignment Operator 
  inline const expdata_ddunit & operator=(const expdata_ddunit & val);
  
  inline expdata_ddunit(const Int_t n, const Double_t ang, 
			const Double_t * e, const Double_t * ee,
			const Double_t * d, const Double_t * ed);
  
  inline expdata_ddunit(const Int_t n, const Double_t ang, 
			const Double_t * e, const Double_t * eel, const Double_t * eeh,
			const Double_t * d, const Double_t * edl, const Double_t * edh);  
  
  inline Double_t GetAngle() const;
  
  inline Double_t GetEnergy(const Int_t i) const;
  inline Double_t GetEnergyError(const Int_t i) const;
  
  inline Double_t GetDifferentialXS(const Int_t i) const;
  inline Double_t GetDifferentialXSError(const Int_t i) const;
  
  TString GetGraphTitle() const;

  inline void ShowYourSelf() const;

  void ShowYourSelf(std::ostringstream & os) const;
  
private:
  Double_t angle;
  
  ClassDef(expdata_ddunit,1)
};

#include "expdata_ddunit.icc"

#endif










