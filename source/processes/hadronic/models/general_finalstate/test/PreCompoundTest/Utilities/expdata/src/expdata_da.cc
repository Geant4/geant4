#include "expdata_da.h"


ClassImp(expdata_da);   

const TGraph * expdata_da::GetGraph() const 
{
  return this->GetData()->GetGraph();
}

void expdata_da::ShowYourSelf(std::ostringstream & os) const 
{
  if (dat) dat->ShowYourSelf(os);
  else os << "There is no data!!\n";
  return;
}

