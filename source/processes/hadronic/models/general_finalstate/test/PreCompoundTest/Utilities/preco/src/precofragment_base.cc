#include "precofragment_base.h"

#include "TClass.h"

ClassImp(precofragment_base)

precofragment_base::precofragment_base() :
  FragmentName("No name"), FragmentA(-1), FragmentZ(-1),
  FragmentM(0.0), FragmentP(0.0,0.0,0.0,0.0),
  FragmentProcess("No process")
{}

precofragment_base::precofragment_base(const precofragment_base& v)
{
  FragmentName = v.FragmentName;
  FragmentA = v.FragmentA;
  FragmentZ = v.FragmentZ;
  FragmentM = v.FragmentM;
  FragmentP = v.FragmentP;
  FragmentProcess = v.FragmentProcess;
}

precofragment_base::precofragment_base(const std::string& frag_name,
				       const int a,
				       const int z,
				       const double m,
				       const TLorentzVector& p,
				       const std::string& process) :
  FragmentName(frag_name), FragmentA(a), FragmentZ(z),
  FragmentM(m), FragmentP(p), FragmentProcess(process)
{}


TBuffer & operator>>(TBuffer & buf, precofragment_base *&obj)
{
  obj = new precofragment_base();
  obj->Streamer(buf);
  delete obj;
  return buf;
}

TBuffer & operator<<(TBuffer & buf, const precofragment_base*  obj)
{
  ((precofragment_base*)obj)->Streamer(buf);
  return buf;
}


