#ifndef G4StatMFParameters_h
#define G4StatMFParameters_h 1

#include "globals.hh"

class G4StatMFParameters
{
private:
  static G4StatMFParameters  theStatMFParameters;

  // +----------------------+
  // | Constant Parameters: |
  // +----------------------+
  // Kappa is used for calculate volume V_f for translational motion of fragments
  static const G4double Kappa;
  // KappaCoulomb is used for calculate Coulomb term energy
  static const G4double KappaCoulomb;
  // Inverse level density
  static const G4double Epsilon0;
  // Bethe-Weizsacker coefficients
  static const G4double E0;
  static const G4double Beta0;
  static const G4double Gamma0;
  // Critical temperature (for liquid-gas phase transitions)
  static const G4double CriticalTemp;
  // Nuclear radius
  static const G4double r0;


  // default constructor
  G4StatMFParameters() 
//    : 
//     Kappa(1.0),
//     KappaCoulomb(2.0),
//     Epsilon0(16.0), // MeV
//     E0(16.0), // MeV
//     Beta0(18.0), // MeV
//     Gamma0(25.0), // MeV
//     CriticalTemp(18.0), // MeV
//     r0(1.17) // fm
    {}


public:

  ~G4StatMFParameters() {};

  static G4StatMFParameters * GetAddress();

  static G4double GetKappa() { return Kappa; }

  static G4double GetKappaCoulomb() { return KappaCoulomb; } 

  static G4double GetEpsilon0() { return Epsilon0; }

  static G4double GetE0() { return E0; }

  static G4double GetBeta0() { return Beta0; } 

  static G4double GetGamma0() { return Gamma0; }

  static G4double GetCriticalTemp() { return CriticalTemp; }

  static G4double Getr0() { return r0; }
};

#endif
