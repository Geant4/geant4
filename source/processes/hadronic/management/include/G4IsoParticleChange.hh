#ifndef G4IsoParticleChange_h
#define G4IsoParticleChange_h

#include "G4Nucleus.hh"

class G4IsoParticleChange 
{
public:
  
  void SetIsotope(const G4String & anIsotope) {theIsotope = anIsotope;}
  void SetProductionPosition(const G4ThreeVector & aPosition) {thePosition = aPosition;}
  void SetProductionTime(const G4double & aProductionTime) {theProductionTime = aProductionTime; }
  void SetParentParticle(const G4DynamicParticle & aProjectile) {theProjectile = aProjectile; }
  void SetMotherNucleus(const G4Nucleus & aTarget) {theTarget = aTarget; }
  void SetProducer(const G4String & aProducer) { theProducer = aProducer; }

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
