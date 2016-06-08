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
//
//
// $Id: G4V3DNucleus.hh,v 1.4.8.1 2001/06/28 19:13:28 gunter Exp $
// GEANT4 tag $Name:  $
//
#ifndef G4V3DNucleus_h
#define G4V3DNucleus_h 1

class G4Nucleon;
#include "G4DynamicParticle.hh"
#include "G4Pair.hh"
#include "Randomize.hh"

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
      
  public:
      G4Pair<G4double, G4double> ChooseImpactXandY(G4double maxImpact);
      G4Pair<G4double, G4double> RefetchImpactXandY(){return theImpactParameter;}

  private:
  
    G4Pair<G4double, G4double> theImpactParameter;

};

inline
G4Pair<G4double, G4double> G4V3DNucleus::
ChooseImpactXandY(G4double maxImpact)
{
  G4double x,y;
  do
  {
    x = 2*G4UniformRand() - 1;
    y = 2*G4UniformRand() - 1;
  }
  while(x*x + y*y > 1);
  G4double impactX = x*(maxImpact); 
  G4double impactY = y*(maxImpact);
  theImpactParameter.first = impactX;
  theImpactParameter.second = impactY;
  return theImpactParameter;
}


#endif


