#include "expdata_unit.h"
#include <stdlib.h>
#include <iostream>

ClassImp(expdata_unit);   


// Copy constructor
expdata_unit::expdata_unit(const expdata_unit & val) : expdata(val)
{
  Xdata = val.Xdata;
  Ydata = val.Ydata;
  Xerror = val.Xerror;
  Yerror = val.Yerror;
  XerrorL = val.XerrorL;
  YerrorL = val.YerrorL;
}

// Symmetric constructor
expdata_unit::expdata_unit(const Int_t n, 
			   const Double_t * x, const Double_t * ex,
			   const Double_t * y, const Double_t * ey)
{
  if (n > 0) 
    {
      Xdata.Set(n,x);
      Ydata.Set(n,y);
      Xerror.Set(n,ex);
      Yerror.Set(n,ey);
      XerrorL.Set(0);
      YerrorL.Set(0);
    }
  else std::cout << "ERROR: Cannot Initialize expdata_unit\n";
}

// Asymmetric constructor
expdata_unit::expdata_unit(const Int_t n, 
			   const Double_t * x, const Double_t * exl, const Double_t * exh,
			   const Double_t * y, const Double_t * eyl, const Double_t * eyh)
{
  if (n > 0) 
    {
      Xdata.Set(n,x);
      Ydata.Set(n,y);
      Xerror.Set(n,exh);
      Yerror.Set(n,eyh);
      XerrorL.Set(n,exl);
      YerrorL.Set(n,eyl);
    }
  else std::cout << "ERROR: Cannot Initialize expdata_unit\n";
}

// Asignment Operator
const expdata_unit & expdata_unit::operator=(const expdata_unit & val)
{
  expdata::operator=(val);
  Xdata = val.Xdata;
  Ydata = val.Ydata;
  Xerror = val.Xerror;
  Yerror = val.Yerror;
  XerrorL = val.XerrorL;
  YerrorL = val.YerrorL;
  return *this;
}


Double_t expdata_unit::GetDataX(const Int_t i) const 
{
  if (!this->IsEmpty()) 
    {
      if (i > 0 && i < this->GetN()) 
	{
	  return Xdata.At(i);
	}    
      else 
	{
	  std::cout << "WARNING: expdata_unit::GetDataX(" << i 
		    << ") out of bounds" << std::endl;

	}
    }
  else 
    {
      std::cout << "WARNING: expdata_unit::GetDataX() has no data!" << std::endl;
    }
  return 0.0;
}



Double_t expdata_unit::GetErrorX(const Int_t i) const
{
  if (!this->IsEmpty()) 
    {
      if (i > 0 && i < this->GetN()) 
	{
	  return Xerror.At(i);
	}    
      else 
	{
	  std::cout << "WARNING: expdata_unit::GetErrorX(" << i 
		    << ") out of bounds" << std::endl;
	    
	}
    }
  else 
    {
      std::cout << "WARNING: expdata_unit::GetErrorX() has no data!" << std::endl;
    }
  return 0.0;
}

Double_t expdata_unit::GetDataY(const Int_t i) const 
{
  if (!this->IsEmpty()) 
    {
      if (i > 0 && i < this->GetN()) 
	{
	  return Ydata.At(i);
	}    
      else 
	{
	  std::cout << "WARNING: expdata_unit::GetDataY(" << i 
		    << ") out of bounds" << std::endl;

	}
    }
  else 
    {
      std::cout << "WARNING: expdata_unit::GetDataY() has no data!" << std::endl;
    }
  return 0.0;
}


Double_t expdata_unit::GetErrorY(const Int_t i) const
{
  if (!this->IsEmpty()) 
    {
      if (i > 0 && i < this->GetN()) 
	{
	  return Yerror.At(i);
	}    
      else 
	{
	  std::cout << "WARNING: expdata_unit::GetErrorY(" << i 
		    << ") out of bounds" << std::endl;

	}
    }
  else 
    {
      std::cout << "WARNING: expdata_unit::GetErrorY() has no data!" << std::endl;
    }
  return 0.0;
}



void expdata_unit::SetData(const Int_t n, 
                           const Double_t * x, const Double_t * ex,
                           const Double_t * y, const Double_t * ey)
{
  this->Reset();
  if (n > 0) 
    {
      Xdata.Set(n,x);
      Ydata.Set(n,y);
      Xerror.Set(n,ex);
      Yerror.Set(n,ey);
      XerrorL.Set(0);
      YerrorL.Set(0);
    }
  return;
}

void expdata_unit::SetData(const Int_t n, const Double_t * x, const Double_t * exl, 
			   const Double_t * exh,
                           const Double_t * y, const Double_t * eyl, 
			   const Double_t * eyh)
{
  this->Reset();
  if (n > 0) 
    {
      Xdata.Set(n,x);
      Ydata.Set(n,y);
      Xerror.Set(n,exh);
      Yerror.Set(n,eyh);
      XerrorL.Set(n,exl);
      YerrorL.Set(n,eyl);
    }
  return;
}

void expdata_unit::Reset() 
{
  Xdata.Set(0);
  Ydata.Set(0);
  Xerror.Set(0);
  Yerror.Set(0);
  XerrorL.Set(0);
  YerrorL.Set(0);
  return;
}


TGraph *  expdata_unit::GetGraph() const
{
  if (this->GetN() <=  0) return 0;
  TGraph * result;
  if (this->IsSymmetric()) 
    {
      result = new TGraphErrors(this->GetN(),Xdata.GetArray(),Ydata.GetArray(),
				Xerror.GetArray(),Yerror.GetArray());
    }
  else 
    {
      result = new TGraphAsymmErrors(this->GetN(),Xdata.GetArray(),Ydata.GetArray(),XerrorL.GetArray(),
				     Xerror.GetArray(),YerrorL.GetArray(),Yerror.GetArray());
    }
  return result;
}





