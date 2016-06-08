// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

#ifndef G4PreCompoundTriton_h
#define G4PreCompoundTriton_h 1

#include "G4VPreCompoundIon.hh"
#include "G4Triton.hh"

class G4PreCompoundTriton : public G4VPreCompoundIon
{
public:
	// default constructor
	G4PreCompoundTriton():G4VPreCompoundIon(3,1) {}

	// copy constructor
	G4PreCompoundTriton(const G4PreCompoundTriton &right): G4VPreCompoundIon(right) {}
	
	// destructor
	~G4PreCompoundTriton() {}

	// operators  
	const G4PreCompoundTriton & operator=(const G4PreCompoundTriton &right) {
		if (&right != this) this->G4VPreCompoundIon::operator=(right);
		return *this;
	}

	G4bool operator==(const G4PreCompoundTriton &right) const
	{ return G4VPreCompoundIon::operator==(right);}

  
	G4bool operator!=(const G4PreCompoundTriton &right) const
	{ return G4VPreCompoundIon::operator!=(right);}

	const G4DynamicParticle GetDynamicParticle() const {
		G4DynamicParticle theDynamicParticle(G4Triton::TritonDefinition(),GetMomentum());
		return theDynamicParticle;
    }




public:
	void CalcExcitonLevelDensityRatios(const G4double Excitons,
												  const G4double Particles)
	{
	// Level density ratios are calculated according to the formula
	// (P!*(N-1)!)/((P-Af)!*(N-1-Af)!*Af! (Af-1)!)
	// where  P is number of particles
	//        N is number of excitons
	//        Af atomic number of emitting fragment
	// the next is a simplification for tritons (Af = 3)

		SetExcitonLevelDensityRatio(((Particles*(Excitons-1.0))*
				 								((Particles-1.0)*(Excitons-2.0)/2.0)*
				 								((Particles-2.0)*(Excitons-3.0)/3.0)/2.0));
   }


	void CalcCondensationProbability(const G4double A)
		// This method computes condensation probability to create a fragment
		// consisting from N nucleons inside a nucleus with A nucleons 
		// This value comes from the formula N^3 (N/A)^(N-1) with N = 3 (triton)
	{
    		SetCondensationProbability(243.0/(A*A));
	}


private:

	virtual G4double GetBarrierPenetrationFactor(const G4double aZ) const;
	virtual G4double GetCCoef(const G4double aZ) const;


};


 
inline G4double G4PreCompoundTriton::GetBarrierPenetrationFactor(const G4double aZ) const
{
	G4double K = 1.0;
	if (aZ>=70.0) {
		K = 0.80;
	} else {
		K = (((0.2357e-5*aZ) - 0.42679e-3)*aZ + 0.27035e-1)*aZ + 0.19025;
	}
	return K+0.12;      
}

inline G4double G4PreCompoundTriton::GetCCoef(const G4double aZ) const
{
	G4double C = 0.0;

	if (aZ >= 70) {
		C = 0.10;
	} else {
		C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
	}

	return C/3.0;                               
}


#endif
