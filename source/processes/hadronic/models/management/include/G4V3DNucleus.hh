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
// $Id: G4V3DNucleus.hh 93818 2015-11-02 11:34:19Z gcosmo $
//
// 20110805  M. Kelsey -- Change nucleon vector to use objects, not pointers

#ifndef G4V3DNucleus_h
#define G4V3DNucleus_h 1

class G4Nucleon;
class G4VNuclearDensity;
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
#include <utility>  
#include <vector>

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
      virtual void Init(G4int theA, G4int theZ) = 0;
      virtual G4bool StartLoop() = 0;
      virtual G4Nucleon * GetNextNucleon() = 0;
      virtual const std::vector<G4Nucleon> & GetNucleons() = 0;
      virtual G4int GetMassNumber() = 0;
      virtual G4double GetMass() = 0;
      virtual G4int GetCharge() = 0;
      virtual G4double GetNuclearRadius() = 0;
      virtual G4double GetNuclearRadius(const G4double maxRelativeDensity) = 0;
      virtual G4double GetOuterRadius() = 0;
      virtual G4double CoulombBarrier() = 0;
      virtual void DoLorentzBoost(const G4LorentzVector & theBoost) = 0;
      virtual void DoLorentzBoost(const G4ThreeVector & theBeta) = 0;
      virtual void DoLorentzContraction(const G4LorentzVector & theBoost) = 0;
      virtual void DoLorentzContraction(const G4ThreeVector & theBeta) = 0;
      virtual void DoTranslation(const G4ThreeVector & theShift) = 0;
      virtual const G4VNuclearDensity * GetNuclearDensity() const = 0;
      virtual void SortNucleonsIncZ() = 0;
      virtual void SortNucleonsDecZ() = 0;

  public:
      std::pair<G4double, G4double> ChooseImpactXandY(G4double maxImpact);
      std::pair<G4double, G4double> RefetchImpactXandY(){return theImpactParameter;}

  private:
  
    std::pair<G4double, G4double> theImpactParameter;

};

inline
std::pair<G4double, G4double> G4V3DNucleus::
ChooseImpactXandY(G4double maxImpact)
{
  G4double x,y;
  do
  {
    x = 2*G4UniformRand() - 1;
    y = 2*G4UniformRand() - 1;
  }
  while(x*x + y*y > 1);  /* Loop checking, 30-Oct-2015, G.Folger */

  G4double impactX = x*(maxImpact); 
  G4double impactY = y*(maxImpact);
  theImpactParameter.first = impactX;
  theImpactParameter.second = impactY;
  return theImpactParameter;
}


#endif


