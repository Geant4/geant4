#ifndef expdata_ddaunit_h
#define expdata_ddaunit_h

#include "expdata_unit.h"
#include "TGraphErrors.h"
#include <iostream>
#include <sstream>

class expdata_ddaunit : public expdata_unit
{
public:
  // Default constructor
  inline expdata_ddaunit();
  // Copy constructor
  inline expdata_ddaunit(const expdata_ddaunit & val);
  // Destructor
  inline ~expdata_ddaunit();
  // Asignement Operator
  inline const expdata_ddaunit & operator=(const expdata_ddaunit & val);

  inline expdata_ddaunit(const Int_t n, const Double_t Emin, const Double_t Emax, 
			 const Double_t * a, const Double_t * ea, 
			 const Double_t * d, const Double_t * ed);
  
  inline expdata_ddaunit(const Int_t n, const Double_t Emin, const Double_t Emax, 
			 const Double_t * a, const Double_t * eal, const Double_t * eah,
			 const Double_t * d, const Double_t * edl, const Double_t * edh);

  inline Double_t GetLowElimit() const;
  inline Double_t GetHighElimit() const;

  inline Double_t GetAngle(const Int_t i) const;
  inline Double_t GetAngleError(const Int_t i) const;

  inline Double_t GetDifferentialXS(const Int_t i) const;
  inline Double_t GetDifferentialXSError(const Int_t i) const;

  TString GetGraphTitle() const;

  void ShowYourSelf(std::ostringstream & os) const;

  inline void ShowYourSelf() const;

private:
  Double_t E_min;
  Double_t E_max;

  ClassDef(expdata_ddaunit,1)
};

#include "expdata_ddaunit.icc"

#endif
