// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara 

#ifndef G4PreCompoundDeuteron_h
#define G4PreCompoundDeuteron_h 1

#include "G4VPreCompoundIon.hh"
#include "G4Deuteron.hh"

class G4PreCompoundDeuteron : public G4VPreCompoundIon
{
public:
	// default constructor
	G4PreCompoundDeuteron():G4VPreCompoundIon(2,1) {};
	// copy constructor
	G4PreCompoundDeuteron(const G4PreCompoundDeuteron &right):
    G4VPreCompoundIon(right) {}
	
	// destructor
	~G4PreCompoundDeuteron() {}

	// operators  
	const G4PreCompoundDeuteron & operator=(const G4PreCompoundDeuteron &right) {
		if (&right != this) this->G4VPreCompoundIon::operator=(right);
		return *this;
	}

	G4bool operator==(const G4PreCompoundDeuteron &right) const
	{ return G4VPreCompoundIon::operator==(right);}
  
	G4bool operator!=(const G4PreCompoundDeuteron &right) const
	{ return G4VPreCompoundIon::operator!=(right);}

	const G4DynamicParticle GetDynamicParticle() const
	{
		G4DynamicParticle theDynamicParticle(G4Deuteron::DeuteronDefinition(),GetMomentum());
		return theDynamicParticle;
	}


public:
	void CalcExcitonLevelDensityRatios(const G4double Excitons,
			             const G4double Particles)
	{
    // Level density ratios are calculated according to the formula
    // (P!*(N-1)!)/((P-Af)!*(N-1-Af)!*Af!)
    // where  P is number of particles
    //        N is number of excitons
    //        Af atomic number of emitting fragment
    // the next is a simplification for deuterons (Af = 2)

		SetExcitonLevelDensityRatio(Particles*(Excitons-1.0)*
				(Particles-1.0)*(Excitons-2.0)/2.0);
   }


	void CalcCondensationProbability(const G4double A)
		// This method computes condensation probability to create a fragment
    // consisting from N nucleons inside a nucleus with A nucleons 
    // This value comes from the formula N^3 (N/A)^(N-1) with N = 2 (deuteron)
	{
		SetCondensationProbability(16.0/A);
	}

private:

	virtual G4double GetBarrierPenetrationFactor(const G4double aZ) const;
	virtual G4double GetCCoef(const G4double aZ) const;

};

inline G4double G4PreCompoundDeuteron::GetBarrierPenetrationFactor(const G4double aZ) const
{
	G4double K = 1.0;
	if (aZ>=70.0) {
		K = 0.80;
	} else {
		K = (((0.2357e-5*aZ) - 0.42679e-3)*aZ + 0.27035e-1)*aZ + 0.19025;
	}
	return K+0.06;                
}

inline G4double G4PreCompoundDeuteron::GetCCoef(const G4double aZ) const
{
	G4double C = 0.0;

	if (aZ >= 70) {
		C = 0.10;
	} else {
		C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
	}
	return C/2.0;
}


#endif
 
