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
#ifndef G4DetailedBalancePhaseSpaceIntegral_h
#define G4DetailedBalancePhaseSpaceIntegral_h

#include "globals.hh"
#include "G4ParticleDefinition.hh"

class G4DetailedBalancePhaseSpaceIntegral
{
  public:
  G4DetailedBalancePhaseSpaceIntegral(const G4ParticleDefinition * aR);
  G4double GetPhaseSpaceIntegral(G4double sqrts);
  
  private:
//  const G4ParticleDefinition * theR;
  const G4double * data;
  
  private:
  static const G4double sqrts[120];
  
  static const G4double delta[120];
  static const G4double delta1600[120];
  static const G4double delta1620[120];
  static const G4double delta1700[120];
  static const G4double delta1900[120];
  static const G4double delta1905[120];
  static const G4double delta1910[120];
  static const G4double delta1920[120];
  static const G4double delta1930[120];
  static const G4double delta1950[120];
  static const G4double N1440[120];
  static const G4double N1520[120];
  static const G4double N1535[120];
  static const G4double N1650[120];
  static const G4double N1675[120];
  static const G4double N1680[120];
  static const G4double N1700[120];
  static const G4double N1710[120];
  static const G4double N1720[120];
  static const G4double N1900[120];
  static const G4double N1990[120];
  static const G4double N2090[120];
  static const G4double N2190[120];
  static const G4double N2220[120];
  static const G4double N2250[120];
};

#endif
