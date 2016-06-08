//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
 













