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
// Calculation of the total, elastic and inelastic cross-sections
// based on parametrisations of (proton, pion, kaon, photon) nucleon
// cross-sections and the hadron-nucleous cross-section model in 
// the framework of Glauber-Gribov approach
//
// 14.03.07 V. Grichine - first implementation
// 04.11.11 V. Grichine - update for kaon-(p,n) xsc, vector spline
// 21.02.12 V. Grichine - update for pion-(p,n) xsc, NS fit++, vector spline
// 30.07.18 V. Ivanchenko - general clean-up
// 30.09.18 V. Grichine hyperon-nucleon xsc first implementation


#ifndef G4HadronNucleonXsc_h
#define G4HadronNucleonXsc_h

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"

class G4Pow;
class G4Element;

class G4HadronNucleonXsc 
{
public:

  G4HadronNucleonXsc ();
  ~G4HadronNucleonXsc ();
   
  // Xsc parametrisations return total x-section
  G4double HadronNucleonXscPDG2005(const G4ParticleDefinition* theParticle, 
			       const G4ParticleDefinition* nucleon, G4double ekin);

  G4double HadronNucleonXscPDG(const G4ParticleDefinition* theParticle, 
			       const G4ParticleDefinition* nucleon, G4double ekin);

  G4double HadronNucleonXscNS(const G4ParticleDefinition* theParticle, 
			      const G4ParticleDefinition* nucleon, G4double ekin);

  G4double KaonNucleonXscNS(const G4ParticleDefinition* theParticle, 
			    const G4ParticleDefinition* nucleon, G4double ekin);

  G4double KaonNucleonXscGG(const G4ParticleDefinition* theParticle, 
			    const G4ParticleDefinition* nucleon, G4double ekin);

  G4double KaonNucleonXscVG(const G4ParticleDefinition* theParticle, 
			    const G4ParticleDefinition* nucleon, G4double ekin);

  G4double HadronNucleonXscVU(const G4ParticleDefinition* theParticle, 
			      const G4ParticleDefinition* nucleon, G4double ekin);

  G4double HadronNucleonXscEL(const G4ParticleDefinition* theParticle, 
			      const G4ParticleDefinition* nucleon, G4double ekin);

  G4double CoulombBarrier(const G4ParticleDefinition* theParticle, 
			  const G4ParticleDefinition* nucleon, G4double ekin);

  // Xsc for G4DynamicParticle projectile
  inline G4double GetHadronNucleonXscEL(const G4DynamicParticle* dp, 
					const G4ParticleDefinition* p)
  { return HadronNucleonXscEL(dp->GetDefinition(), p, dp->GetKineticEnergy()); }

  inline G4double GetHadronNucleonXscPDG(const G4DynamicParticle* dp, 
					 const G4ParticleDefinition* p)
  { return HadronNucleonXscPDG(dp->GetDefinition(), p, dp->GetKineticEnergy()); }

  inline G4double GetHadronNucleonXscNS(const G4DynamicParticle* dp, 
					const G4ParticleDefinition* p)
  { return HadronNucleonXscNS(dp->GetDefinition(), p, dp->GetKineticEnergy()); }

  inline G4double GetKaonNucleonXscGG(const G4DynamicParticle* dp, 
				      const G4ParticleDefinition* p)
  { return KaonNucleonXscGG(dp->GetDefinition(), p, dp->GetKineticEnergy()); }

  inline G4double GetHadronNucleonXscVU(const G4DynamicParticle* dp, 
					const G4ParticleDefinition* p)
  { return HadronNucleonXscVU(dp->GetDefinition(), p, dp->GetKineticEnergy()); }

  inline G4double GetCoulombBarrier(const G4DynamicParticle* dp, 
				    const G4ParticleDefinition* p)
  { return CoulombBarrier(dp->GetDefinition(), p, dp->GetKineticEnergy()); }

  // Xsc access
  inline G4double GetTotalHadronNucleonXsc()     const { return fTotalXsc;     }; 
  inline G4double GetElasticHadronNucleonXsc()   const { return fElasticXsc;   }; 
  inline G4double GetInelasticHadronNucleonXsc() const { return fInelasticXsc; }; 

  // obsolete methods
  G4bool IsApplicable(const G4DynamicParticle* aDP, const G4Element*);

  G4bool IsIsoApplicable(const G4DynamicParticle* aDP, G4int Z);

  void DumpPhysicsTable(const G4ParticleDefinition&) 
  {G4cout << "G4HadronNucleonXsc: uses parametrisation"<<G4endl;}

  void CrossSectionDescription(std::ostream&) const;

private:

  G4double ComputeKaonNucleonXsc(const G4ParticleDefinition* theParticle, 
				 const G4ParticleDefinition* nucleon, G4double ekin);


  inline G4double CalcMandelstamS(G4double ekin1, G4double mass1, G4double mass2)
  { return mass1*mass1 + mass2*mass2 + 2*mass2*(ekin1 + mass1); }

  inline G4double CalculateEcmValue(G4double ekin1, G4double mass1, G4double mass2)
  { return std::sqrt(CalcMandelstamS(ekin1, mass1, mass2)); };

  G4double fLowerLimit; 
  G4double fTotalXsc, fElasticXsc, fInelasticXsc;
  G4double fHypTotXscCof;
  G4Pow* g4calc;

  const G4ParticleDefinition* theGamma;
  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theNeutron;
  const G4ParticleDefinition* theAProton;
  const G4ParticleDefinition* theANeutron;
  const G4ParticleDefinition* thePiPlus;
  const G4ParticleDefinition* thePiMinus;
  const G4ParticleDefinition* thePiZero;
  const G4ParticleDefinition* theKPlus;
  const G4ParticleDefinition* theKMinus;
  const G4ParticleDefinition* theK0S;
  const G4ParticleDefinition* theK0L;
  const G4ParticleDefinition* theL;
  const G4ParticleDefinition* theAntiL;
  const G4ParticleDefinition* theSPlus;
  const G4ParticleDefinition* theASPlus;
  const G4ParticleDefinition* theSMinus;
  const G4ParticleDefinition* theASMinus;
  const G4ParticleDefinition* theS0;
  const G4ParticleDefinition* theAS0;
  const G4ParticleDefinition* theXiMinus;
  const G4ParticleDefinition* theXi0;
  const G4ParticleDefinition* theAXiMinus;
  const G4ParticleDefinition* theAXi0;
  const G4ParticleDefinition* theOmega;
  const G4ParticleDefinition* theAOmega;
  const G4ParticleDefinition* theD;
  const G4ParticleDefinition* theT;
  const G4ParticleDefinition* theA;
  const G4ParticleDefinition* theHe3;
};


#endif
