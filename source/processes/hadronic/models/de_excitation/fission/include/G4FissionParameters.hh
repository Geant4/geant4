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
//
// $Id: G4FissionParameters.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4FissionParameters_h
#define G4FissionParameters_h 1

#include "globals.hh"


class G4FissionParameters 
{
public:
  // Only available constructor
  G4FissionParameters(G4int A, G4int Z, G4double ExEnergy, G4double FissionBarrier);

  ~G4FissionParameters();  

private:  
  // Default constructor
  G4FissionParameters();

  // Copy constructor
  G4FissionParameters(const G4FissionParameters &right);

  const G4FissionParameters & operator=(const G4FissionParameters &right);
  G4bool operator==(const G4FissionParameters &right) const;
  G4bool operator!=(const G4FissionParameters &right) const;
  
public:

  inline G4double GetA1(void) const { return A1; }
  inline G4double GetA2(void) const { return A2; }

  inline G4double GetAs(void) const { return As; }
  inline G4double GetSigma1(void) const { return Sigma1; }
  inline G4double GetSigma2(void) const { return Sigma2; }
  inline G4double GetSigmaS(void) const { return SigmaS; }
  inline G4double GetW(void) const { return w; }

private:

  // Mean numbers of the corresponding Gaussians for assymmetric
  // fission
  static const G4double A1;
  static const G4double A2;

  // Mean number for symmetric fission
  G4double As;

  // Dispersions of the corresponding Gaussians for assymmetric
  // fission
  G4double Sigma1;
  G4double Sigma2;

  // Dispersion for symmetric fission
  G4double SigmaS;

  // Weight which determines the relative contribution of symmetric
  // and assymmetric components
  G4double w;


};


#endif
