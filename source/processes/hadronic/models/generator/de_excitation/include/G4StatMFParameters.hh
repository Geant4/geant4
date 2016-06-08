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
  static const G4double _Kappa;
  // KappaCoulomb is used for calculate Coulomb term energy
  static const G4double _KappaCoulomb;
  // Inverse level density
  static const G4double _Epsilon0;
  // Bethe-Weizsacker coefficients
  static const G4double _E0;
  static const G4double _Beta0;
  static const G4double _Gamma0;
  // Critical temperature (for liquid-gas phase transitions)
  static const G4double _CriticalTemp;
  // Nuclear radius
  static const G4double _r0;


  // default constructor
  G4StatMFParameters() 
//    : 
//     _Kappa(1.0),
//     _KappaCoulomb(2.0),
//     _Epsilon0(16.0), // MeV
//     _E0(16.0), // MeV
//     _Beta0(18.0), // MeV
//     _Gamma0(25.0), // MeV
//     _CriticalTemp(18.0), // MeV
//     _r0(1.17) // fm
    {}


public:

  ~G4StatMFParameters() {};

  static G4StatMFParameters * GetAddress();

  static G4double GetKappa() { return _Kappa; }

  static G4double GetKappaCoulomb() { return _KappaCoulomb; } 

  static G4double GetEpsilon0() { return _Epsilon0; }

  static G4double GetE0() { return _E0; }

  static G4double GetBeta0() { return _Beta0; } 

  static G4double GetGamma0() { return _Gamma0; }

  static G4double GetCriticalTemp() { return _CriticalTemp; }

  static G4double Getr0() { return _r0; }
  
  static G4double Beta(const G4double T);
  
  static G4double DBetaDT(const G4double T);
};

#endif
