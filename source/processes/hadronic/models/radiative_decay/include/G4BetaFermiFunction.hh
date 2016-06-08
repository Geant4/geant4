#ifndef G4BetaFermiFunction_h
#define G4BetaFermiFunction_h 1

#include "globals.hh"

class G4BetaFermiFunction
{
  // class  description
  // It is to calculate the Coulomb correction to beta particles
  // 

public: // with description
  
  G4BetaFermiFunction(G4int const fA, G4int const fZ) :
    A(fA), Z(fZ)
  {};
  // constructor: fA the daughter nucleus mass.
  //              fZ the daughter nucleus charge. Negative value for 
  //                 beta+ decays.
  //
  ~G4BetaFermiFunction() 
  // desctructor
  //
  {};
  G4double GetFF(const G4double E);
  // Returns the BetaFermi factor at energy E.
  // E kinetic energy of the beta particle in unit of Me.
  //
  G4double GetFFN(const G4double E0);
  // Returns the BetaFermi factor normalisation, i.e. the maximum
  // value for beta decay with end-point energy E0.
  // E0 dose not including beta particle rest mass in unit
  // of Me.
  //
private:

  G4int A;
  G4int Z;
  
  static const G4double PI;

private:

  G4double Gamma(G4double X);
  
};
#endif
 













