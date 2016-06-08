// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara 

#ifndef G4PreCompoundAlpha_h
#define G4PreCompoundAlpha_h 1

#include "G4VPreCompoundIon.hh"
#include "G4Alpha.hh"


class G4PreCompoundAlpha : public G4VPreCompoundIon
{
public:
	// default constructor
	G4PreCompoundAlpha():G4VPreCompoundIon(4,2) {}

	// copy constructor
	G4PreCompoundAlpha(const G4PreCompoundAlpha &right): G4VPreCompoundIon(right) {}

	// destructor
	~G4PreCompoundAlpha() {}

	// operators  
	const G4PreCompoundAlpha & operator=(const G4PreCompoundAlpha &right) {
		if (&right != this) this->G4VPreCompoundIon::operator=(right);
		return *this;
	}

	G4bool operator==(const G4PreCompoundAlpha &right) const
	{ return G4VPreCompoundIon::operator==(right);}
  
	G4bool operator!=(const G4PreCompoundAlpha &right) const
	{ return G4VPreCompoundIon::operator!=(right);}

	const G4DynamicParticle GetDynamicParticle() const {
		G4DynamicParticle theDynamicParticle(G4Alpha::AlphaDefinition(),GetMomentum());
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
	// the next is a simplification for alphas (Af = 4)
    
		SetExcitonLevelDensityRatio(((Particles*(Excitons-1.0))*
											((Particles-1.0)*(Excitons-2.0)/2.0)*
				 							((Particles-2.0)*(Excitons-3.0)/3.0)*
						 					((Particles-3.0)*(Excitons-4.0)/4.0))/6.0);
	}



	void CalcCondensationProbability(const G4double A)
	// This method computes condensation probability to create a fragment
	// consisting from N nucleons inside a nucleus with A nucleons 
	// This value comes from the formula N^3 (N/A)^(N-1) with N = 4 (alpha)
	{
		SetCondensationProbability(4096.0/(A*A*A));
	}

private:

	virtual G4double GetBarrierPenetrationFactor(const G4double aZ) const;
	virtual G4double GetCCoef(const G4double aZ) const;

};

inline G4double G4PreCompoundAlpha::GetBarrierPenetrationFactor(const G4double aZ) const
{
	G4double K = 1.0;
	if (aZ>=70.0) {
		K = 0.98;
	} else {
		K = (((0.23684e-5*aZ) - 0.42143e-3)*aZ + 0.25222e-1)*aZ + 0.46699;
	}
	return K;    
}

inline G4double G4PreCompoundAlpha::GetCCoef(const G4double aZ) const
{
	G4double C = 0.0;

	if (aZ <= 30) {
		C = 0.10;
	} else if (aZ <= 50) {
		C = 0.1 + -((aZ-50.)/20.)*0.02;
	} else if (aZ < 70) {
		C = 0.08 + -((aZ-70.)/20.)*0.02;
	} else {
		C = 0.06;
	}
	return C;            
}

#endif
 
