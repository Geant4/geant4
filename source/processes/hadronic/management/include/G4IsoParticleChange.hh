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
#ifndef G4IsoParticleChange_h
#define G4IsoParticleChange_h

#include "G4Nucleus.hh"
// Class Description
// THis is the class you can retrieve from the hadronic process, that contains
// the isotope production information/
// Class Description - End

class G4IsoParticleChange 
{
public:
  
  void SetIsotope(const G4String & anIsotope) {theIsotope = anIsotope;}
  void SetProductionPosition(const G4ThreeVector & aPosition) {thePosition = aPosition;}
  void SetProductionTime(const G4double & aProductionTime) {theProductionTime = aProductionTime; }
  void SetParentParticle(const G4DynamicParticle & aProjectile) {theProjectile = aProjectile; }
  void SetMotherNucleus(const G4Nucleus & aTarget) {theTarget = aTarget; }
  void SetProducer(const G4String & aProducer) { theProducer = aProducer; }

public:// With description
  // This is the information you can retrieve.
  
  G4String GetIsotope() {return theIsotope;}
  G4ThreeVector GetProductionPosition() {return thePosition;}
  G4double GetProductionTime() {return theProductionTime;}
  G4DynamicParticle GetParentParticle() {return theProjectile;}
  G4Nucleus GetMotherNucleus() {return theTarget;}
  G4String GetProducer() {return theProducer;}

private:

  G4String theIsotope;
  G4ThreeVector thePosition;
  G4double theProductionTime;
  G4DynamicParticle theProjectile;
  G4Nucleus theTarget;
  G4String theProducer;
};

#endif
