#ifndef precofragment_h
#define precofragment_h

#include "TObject.h"
#include "precofragment_base.h"

class precofragment : public TObject, public precofragment_base
{
public:
  precofragment() : precofragment_base() {}
  precofragment(const std::string& frag_name,
		const int a,
		const int z,
		const double m,
		const TLorentzVector& p,
		const std::string& process) :
    precofragment_base(frag_name,a,z,m,p,process) {}

  virtual ~precofragment() {}

ClassDef(precofragment,1)

};
#endif
