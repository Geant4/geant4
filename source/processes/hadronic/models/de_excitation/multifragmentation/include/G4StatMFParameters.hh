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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

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
  
  static G4double GetMaxAverageMultiplicity(const G4int A);
};

#endif
