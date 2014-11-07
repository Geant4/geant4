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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPKallbachMannSyst_h
#define G4ParticleHPKallbachMannSyst_h 1

#include "globals.hh"

class G4ParticleHPKallbachMannSyst
{
  public:
  
   G4ParticleHPKallbachMannSyst(G4double aCompoundFraction,
                               G4double anIncidentEnergy, G4double anIncidentMass,
                               G4double aProductEnergy, G4double aProductMass,
                               G4double aResidualMass, G4int aResidualA, G4int aResidualZ,
                               G4double aTargetMass, G4int aTargetA, G4int aTargetZ)
  {
    theCompoundFraction = aCompoundFraction;
    theIncidentEnergy = anIncidentEnergy;
    theIncidentMass = anIncidentMass;
    theProductEnergy = aProductEnergy;
    theProductMass = aProductMass;
    theResidualMass = aResidualMass;
    theResidualA = aResidualA;
    theResidualZ = aResidualZ;
    theTargetMass = aTargetMass;
    theTargetA = aTargetA;
    theTargetZ = aTargetZ;
  }
  
  ~G4ParticleHPKallbachMannSyst() {};
  
  G4double Sample(G4double anEnergy);
  
  G4double Kallbach(G4double cosTh, G4double anEnergy);
  
  G4double GetKallbachZero(G4double anEnergy);
  
  G4double A(G4double anEnergy);
  
  G4double SeparationEnergy(G4int Ac, G4int Nc, G4int AA, G4int ZA);
  
  private:
  
  G4double theCompoundFraction;
  G4double theIncidentEnergy;
  G4double theIncidentMass;
  G4double theProductEnergy;
  G4double theProductMass;
  G4double theResidualMass;
  G4double theTargetMass;
  G4int theResidualA;
  G4int theResidualZ;
  G4int theTargetA;
  G4int theTargetZ;
};

#endif
