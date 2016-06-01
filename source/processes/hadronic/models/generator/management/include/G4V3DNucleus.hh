// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4V3DNucleus.hh,v 1.2 1998/10/30 16:36:36 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
#ifndef G4V3DNucleus_h
#define G4V3DNucleus_h 1

#include "G4Nucleon.hh"
#include "G4DynamicParticle.hh"

class G4V3DNucleus 
{

  public:
      G4V3DNucleus();
      virtual ~G4V3DNucleus();

  private:
      G4V3DNucleus(const G4V3DNucleus &right);
      const G4V3DNucleus & operator=(const G4V3DNucleus &right);
      int operator==(const G4V3DNucleus &right) const;
      int operator!=(const G4V3DNucleus &right) const;

  public:
      virtual void Init(G4double theA, G4double theZ) = 0;
      virtual G4bool StartLoop() = 0;
      virtual G4Nucleon * GetNextNucleon() = 0;
      virtual G4int GetMassNumber() = 0;
      virtual G4double GetMass() = 0;
      virtual G4int GetCharge() = 0;
      virtual G4double GetNuclearRadius() = 0;
      virtual G4double GetNuclearRadius(const G4double maxRelativeDensity) = 0;
      virtual G4double GetOuterRadius() = 0;
      virtual void DoLorentzBoost(const G4LorentzVector & theBoost) = 0;
      virtual void DoLorentzBoost(const G4ThreeVector & theBeta) = 0;
      virtual void DoLorentzContraction(const G4LorentzVector & theBoost) = 0;
      virtual void DoLorentzContraction(const G4ThreeVector & theBeta) = 0;
      virtual void DoTranslation(const G4ThreeVector & theShift) = 0;

  private:

};

#endif


