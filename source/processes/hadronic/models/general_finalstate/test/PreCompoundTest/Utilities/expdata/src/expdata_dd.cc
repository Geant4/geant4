#include "expdata_dd.h"
#include <stdlib.h>
#include <iomanip>

ClassImp(expdata_dd);

expdata_dd::expdata_dd() : nangles(0)
{
  dat = new TObjArray(5);
}


const TGraph * expdata_dd::GetGraph(const Int_t i) const 
{
  if (i >= 0 && i < nangles) 
    {
      return this->GetData(i)->GetGraph();
    }
  else 
    {
      std::cout << "WARNING!! expdata_dd::GetGraph("<<i<<"): Out of bounds\n";
    }
  return 0;
}

expdata_ddunit * expdata_dd::GetData(const Double_t ang) const
{
  for (Int_t i = 0; i < this->GetNangles(); i++) 
    {
      if (this->GetData(i)->GetAngle() == ang) 
	return this->GetData(i);
    }
  std::cout << "expdata_dd::GetData(" << ang << "): Angle not found\n";
  return 0;
}


void expdata_dd::ShowYourSelf(const Int_t v, std::ostringstream & os) const 
{
  os << "****************************************************************\n"
     << "*               ENERGY-ANGLE DISTRIBUTION                      *\n" 
     << "****************************************************************\n\n";

  if (v < 0 || v > nangles-1) 
    {
      os << "+-------------+\n"
	 << "| " << std::setw(4) << nangles << " Angles |  ";
      for (Int_t i = 0; i < nangles; i++) 
	{
	  os << ((expdata_ddunit*)(this->GetData(i)))->GetAngle() << " ";
	}
      os << "degrees\n"
	 << "+-------------+\n\n";
    }
  else 
    {
      ((expdata_ddunit*)(this->GetData(v)))->ShowYourSelf(os);
    } 
  
  return;
}





