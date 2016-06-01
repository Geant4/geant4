#ifndef G4StatMFFragment_h
#define G4StatMFFragment_h 1


#include "G4StatMFParameters.hh"

class G4StatMFFragment {
public:
  // default constructor
  G4StatMFFragment():
    InvLevelDensity(0.0),
    ZARatio(0.0),
    DegeneracyFactor(0.0),
    Multiplicity(0.0),
    A(0.0),
    Z(0.0),
    Energy(0.0)
    {};
  // destructor
  ~G4StatMFFragment() {};

private:
  // copy constructor
  G4StatMFFragment(const G4StatMFFragment & right);

  // operators
  const G4StatMFFragment & operator=(const G4StatMFFragment & right);

public:
  G4bool operator==(const G4StatMFFragment & right) const;
  G4bool operator!=(const G4StatMFFragment & right) const;



private:
  
  // Inverse Level Density
  G4double InvLevelDensity;

  // Z/A ratio
  G4double ZARatio;


  // Degeneracy Factor
  G4double DegeneracyFactor;

  // Fragments Multiplicitie
  G4double Multiplicity;

  // Atomic number
  G4double A;
  
  // Charge
  G4double Z;


  // Energy
  G4double Energy;

public:

  void SetInvLevelDensity(const G4double value) {
    InvLevelDensity = value;
  }
  void SetInvLevelDensity(const G4int value) {
    // 
    if (value == 0) InvLevelDensity = 0.0;
    else InvLevelDensity = G4StatMFParameters::GetEpsilon0()/
	   (1.0+0.002*((value+1.0)/25.0)*((value+1.0)/25.0));
  }
  const G4double GetInvLevelDensity() const {
    return InvLevelDensity;
  }


  void SetZARatio(const G4double value) {
     ZARatio = value;
  }
  const G4double GetZARatio() const {
    return ZARatio;
  }


  void SetDegeneracyFactor(const G4double value) {
    DegeneracyFactor = value;
  }
  const G4double GetDegeneracyFactor() const {
    return DegeneracyFactor;
  }


  void SetMultiplicity(const G4double value) {
     Multiplicity = value;
  }
  const G4double GetMultiplicity() const {
    return Multiplicity;
  }


  void SetA(const G4double value) {
    A = value;
  }
  const G4double GetA() const {
    return A;
  }

  void SetZ(const G4double value) {
    Z = value;
  }
  const G4double GetZ() const {
    return Z;
  }

  void SetEnergy(const G4double value) {
    Energy = value;
  }
  const G4double GetEnergy() const {
    return Energy;
  }



};



#endif
