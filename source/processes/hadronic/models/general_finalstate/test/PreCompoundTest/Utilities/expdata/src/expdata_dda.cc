#include "expdata_dda.h"
#include <stdlib.h>
#include <iomanip>

ClassImp(expdata_dda);

expdata_dda::expdata_dda() : nranges(0)
{
  dat = new TObjArray(5);
}

const TGraph * expdata_dda::GetGraph(const Int_t i) const
{
  if (i >= 0 && i < nranges)
    {
      return this->GetData(i)->GetGraph();
    }
  else
    {
      std::cout << "WARNING!! expdata_dda::GetGraph("<<i<<"): Out of bounds\n";
    }
  return 0;
}

expdata_ddaunit * expdata_dda::GetData(const Double_t emin, const Double_t emax)
{
  for (Int_t i = 0; i < this->GetNranges(); i++)
    {
      if (this->GetData(i)->GetLowElimit() == emin && this->GetData(i)->GetHighElimit() == emax)
	return  this->GetData(i);
    }
  std::cout << "expdata_dda::GetData("<< emin << ',' << emax <<"): Range not found\n";
  return 0;
}

void expdata_dda::ShowYourSelf(const Int_t v, std::ostringstream & os) const
{
  os << "****************************************************************\n"
     << "*               ENERGY-ANGLE DISTRIBUTION                      *\n" 
     << "****************************************************************\n\n";

  if (v < 0 || v > nranges-1) 
    {
      os << "+-------------+\n"
	 << "| " << std::setw(4) << nranges << " Ranges |  ";
      for (Int_t i = 0; i < nranges; i++) 
	{
	  os << '[' << ((expdata_ddaunit*)(this->GetData(i)))->GetLowElimit()
	     << ',' << ((expdata_ddaunit*)(this->GetData(i)))->GetHighElimit() << "] ";
	}
      os << "MeV\n"
	 << "+-------------+\n\n";
    }
  else 
    {
      ((expdata_ddaunit*)(this->GetData(v)))->ShowYourSelf(os);
    } 
  
  return;
}
