//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
 













