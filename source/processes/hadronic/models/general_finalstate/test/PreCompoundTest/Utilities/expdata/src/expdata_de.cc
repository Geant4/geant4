#include "expdata_de.h"


ClassImp(expdata_de);   

const TGraph * expdata_de::GetGraph() const 
{
  return this->GetData()->GetGraph();
}

void expdata_de::ShowYourSelf(std::ostringstream & os) const 
{
  if (dat) dat->ShowYourSelf(os);
  else os << "There is no data!!\n";
  return;
}

