#ifndef precofragment_base_h
#define precofragment_base_h

#include <string>
#include "TVector3.h"
#include "TLorentzVector.h"

class precofragment_base
{
private:
    std::string    FragmentName;
    int            FragmentA;
    int            FragmentZ;
    double         FragmentM;
    TLorentzVector FragmentP;
    std::string    FragmentProcess;

public:
    precofragment_base();
    precofragment_base(const precofragment_base&);
    precofragment_base(const std::string& frag_name,
		       const int a,
		       const int z,
		       const double m,
		       const TLorentzVector& p,
		       const std::string& process);
  virtual ~precofragment_base() {}
    
    // Get Methods
    inline std::string GetFragmentName() const { return FragmentName; }
    inline std::string GetProcessName() const { return FragmentProcess; }
    inline int GetA() const { return FragmentA; }
    inline int GetZ() const { return FragmentZ; }
    inline double GetMass() const { return FragmentM; }
    inline double GetKineticE() const { return FragmentP.E()-FragmentM; }
    inline double GetKineticE(const TVector3 boost) const 
	{ 
	    TLorentzVector P(FragmentP);
	    P.Boost(boost);
	    return P.E()-FragmentM; 
	}
    inline double GetTotalE() const { return FragmentP.E(); }
    inline TLorentzVector GetMomentum() const { return FragmentP; }
    
    
    inline double GetTheta() const { return FragmentP.Theta(); }
    inline double GetTheta(const TVector3& dir) const
	{
	    return FragmentP.Angle(dir);
	}
    inline double GetTheta(const TVector3& dir, const TVector3& boost) const
	{
	  TLorentzVector P(FragmentP);
	  P.Boost(boost);
	  return P.Angle(dir);
	}
    friend TBuffer &operator<<(TBuffer& b, const precofragment_base * frag);
    friend TBuffer &operator>>(TBuffer& b, const precofragment_base * frag);
    
    ClassDef(precofragment_base,1)

};
#endif
