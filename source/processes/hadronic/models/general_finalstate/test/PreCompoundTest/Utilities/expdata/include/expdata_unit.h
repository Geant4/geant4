#ifndef expdata_unit_h
#define expdata_unit_h

#include "expdata.h"

#include "TObject.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TArrayD.h"
#include <iostream>
#include <sstream>

class expdata_unit : public expdata
{
public:
    
  void Reset();
  
  inline Bool_t IsEmpty() const;
  inline Bool_t IsAsymmetric() const;
  inline Bool_t IsSymmetric() const;
  
  inline Int_t GetN() const;

  void SetData(const Int_t n, const Double_t * X, const Double_t * EX,
	       const Double_t * Y, const Double_t * EY);
  
  void SetData(const Int_t n, const Double_t * X, const Double_t * EXl, const Double_t * EXh,
	       const Double_t * Y, const Double_t * EYl, const Double_t * EYh);
  
  
  virtual TString GetGraphTitle() const = 0;
  
  inline void ShowYourSelf() const;
  inline void ShowYourSelf(std::ostringstream & os) const;
  
  TGraph *  GetGraph() const;
    
protected:
  // Default constructor
  inline expdata_unit();
    
  // Destructor
  inline ~expdata_unit();
    
  // Copy constuctor
  expdata_unit(const expdata_unit & val);
    
  // Symmetric constructor
  expdata_unit(const Int_t n, 
	       const Double_t * x, const Double_t * ex,
	       const Double_t * y, const Double_t * ey);
    
  // Asymmetric constructor
  expdata_unit(const Int_t n, 
	       const Double_t * x, const Double_t * exl, const Double_t * exh,
	       const Double_t * y, const Double_t * eyl, const Double_t * eyh);
    
  // Assignment Operator
  const expdata_unit & operator=(const expdata_unit &);
    
    
  Double_t GetDataX(const Int_t i) const;   
  Double_t GetErrorX(const Int_t i) const;
    
    
  Double_t GetDataY(const Int_t i) const;
  Double_t GetErrorY(const Int_t i) const;
    
private:
   
  TArrayD  Xdata; 
  TArrayD  Ydata; 
  TArrayD  Xerror;
  TArrayD  Yerror;
  TArrayD  XerrorL;
  TArrayD  YerrorL;

  ClassDef(expdata_unit,1)
};

#include "expdata_unit.icc"

#endif
