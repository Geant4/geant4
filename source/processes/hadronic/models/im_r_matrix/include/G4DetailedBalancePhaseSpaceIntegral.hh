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
#ifndef G4DetailedBalancePhaseSpaceIntegral_h
#define G4DetailedBalancePhaseSpaceIntegral_h

#include "globals.hh"
#include "G4ParticleDefinition.hh"

class G4DetailedBalancePhaseSpaceIntegral
{
  public:
  G4DetailedBalancePhaseSpaceIntegral(G4ParticleDefinition * aR);
  G4double GetPhaseSpaceIntegral(G4double sqrts);
  
  private:
  G4ParticleDefinition * theR;
  G4double * data;
  
  private:
  static G4double sqrts[120];
  
  static G4double delta[120];
  static G4double delta1600[120];
  static G4double delta1620[120];
  static G4double delta1700[120];
  static G4double delta1900[120];
  static G4double delta1905[120];
  static G4double delta1910[120];
  static G4double delta1920[120];
  static G4double delta1930[120];
  static G4double delta1950[120];
  static G4double N1440[120];
  static G4double N1520[120];
  static G4double N1535[120];
  static G4double N1650[120];
  static G4double N1675[120];
  static G4double N1680[120];
  static G4double N1700[120];
  static G4double N1710[120];
  static G4double N1720[120];
  static G4double N1900[120];
  static G4double N1990[120];
  static G4double N2090[120];
  static G4double N2190[120];
  static G4double N2220[120];
  static G4double N2250[120];
};

#endif
