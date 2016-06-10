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
// $Id: G4StatMFParameters.hh 91834 2015-08-07 07:24:22Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFParameters_h
#define G4StatMFParameters_h 1

#include "globals.hh"

class G4StatMFParameters
{
public:
  
  G4StatMFParameters();

  ~G4StatMFParameters();
  
  static G4double GetKappa();
  
  static G4double GetKappaCoulomb(); 
  
  static G4double GetEpsilon0();
  
  static G4double GetE0();
  
  static G4double GetBeta0(); 
  
  static G4double GetGamma0();
  
  static G4double GetCriticalTemp();
  
  static G4double Getr0();

  static G4double GetCoulomb();
  
  static G4double Beta(G4double T);
  
  static G4double DBetaDT(G4double T);
  
  static G4double GetMaxAverageMultiplicity(G4int A);

  // +----------------------+
  // | Constant Parameters: |
  // +----------------------+
  // Kappa is used for calculate volume V_f for translational 
  // motion of fragments
  static const G4double fKappa;
  // KappaCoulomb is used for calculate Coulomb term energy
  static const G4double fKappaCoulomb;
  // Inverse level density
  static const G4double fEpsilon0;
  // Bethe-Weizsacker coefficients
  static const G4double fE0;
  static const G4double fBeta0;
  static const G4double fGamma0;
  // Critical temperature (for liquid-gas phase transitions)
  static const G4double fCriticalTemp;
  // Nuclear radius
  static const G4double fr0;
  // Coulomb 
  static const G4double fCoulomb;

};

#endif
