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
// Calculation of the nucleus-nucleus total, inelastic, production, 
// elastic and quasi-elastic  cross-sections
// based on parametrisations of nucleon-nucleon
// cross-sections  in 
// the framework of simplified Glauber-Gribov approach
//
//
// 24.11.08 V. Grichine - first implementation based on 
//                        G4GlauberGribovCrossSection
//
// 04.09.18 V. Ivantchenko Major revision of interfaces and implementation
// 27.05.19 V. Ivantchenko Removed obsolete methods and members 
//

#ifndef G4ComponentGGNuclNuclXsc_h
#define G4ComponentGGNuclNuclXsc_h

#include "globals.hh"
#include "G4VComponentCrossSection.hh"
#include "G4DynamicParticle.hh"

class G4ParticleDefinition;
class G4HadronNucleonXsc;
class G4ComponentGGHadronNucleusXsc;
class G4Material;

class G4ComponentGGNuclNuclXsc : public G4VComponentCrossSection
{
public:

  G4ComponentGGNuclNuclXsc ();
  virtual ~G4ComponentGGNuclNuclXsc ();

  // virtual interface methods
  G4double GetTotalElementCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy, 
				       G4int Z, G4double A) final;

  G4double GetTotalIsotopeCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy,
				       G4int Z, G4int A) final;

  G4double GetInelasticElementCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4double A) final;

  G4double GetInelasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4int A) final;

  G4double GetElasticElementCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4double A) final;

  G4double GetElasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int A) final;
 
  G4double ComputeQuasiElasticRatio(const G4ParticleDefinition* aParticle,
				    G4double kinEnergy, 
				    G4int Z, G4int A) final;

  void BuildPhysicsTable(const G4ParticleDefinition&) final;

  void DumpPhysicsTable(const G4ParticleDefinition&) final;

  void Description(std::ostream&) const final;
   
  // Extra methods
  //  inline G4double GetElementCrossSection(const G4DynamicParticle*, 
  //				         G4int Z, const G4Material*);

  inline G4double GetZandACrossSection(const G4DynamicParticle*, 
				       G4int Z, G4int A);

  inline G4double GetCoulombBarier(const G4DynamicParticle*, 
			           G4double Z, G4double A, 
                                   G4double pR, G4double tR);

  G4double ComputeCoulombBarier(const G4ParticleDefinition* aParticle,
				G4double kinEnergy, G4int Z, G4int A,
				G4double pR, G4double tR);

  G4double GetRatioSD(const G4DynamicParticle*, G4double At, G4double Zt);
  G4double GetRatioQE(const G4DynamicParticle*, G4double At, G4double Zt);

  inline G4double GetElasticGlauberGribov(const G4DynamicParticle*,G4int Z, G4int A);
  inline G4double GetInelasticGlauberGribov(const G4DynamicParticle*,G4int Z, G4int A);

  inline G4double GetTotalGlauberGribovXsc() const       { return fTotalXsc;     }; 
  inline G4double GetElasticGlauberGribovXsc() const     { return fElasticXsc;   }; 
  inline G4double GetInelasticGlauberGribovXsc() const   { return fInelasticXsc; }; 
  inline G4double GetProductionGlauberGribovXsc() const  { return fProductionXsc; }; 
  inline G4double GetDiffractionGlauberGribovXsc() const { return fDiffractionXsc; }; 

private:

  // Glauber-Gribov cross section
  void ComputeCrossSections(const G4ParticleDefinition* aParticle,
			    G4double kinEnergy, G4int Z, G4int A);

  G4double fTotalXsc, fElasticXsc, fInelasticXsc;
  G4double fProductionXsc, fDiffractionXsc;
  // Cache
  G4double fEnergy;
 
  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theNeutron;
  const G4ParticleDefinition* theLambda;

  G4ComponentGGHadronNucleusXsc* fHadrNucl; 
  G4HadronNucleonXsc* fHNXsc;

  // Cache
  const G4ParticleDefinition* fParticle;
  G4int fZ, fA;    
};

inline G4double
G4ComponentGGNuclNuclXsc::GetElasticGlauberGribov(const G4DynamicParticle* dp,
                                                  G4int Z, G4int A)
{
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(), Z, A);
  return fElasticXsc;
}

inline G4double
G4ComponentGGNuclNuclXsc::GetInelasticGlauberGribov(const G4DynamicParticle* dp,
                                                    G4int Z, G4int A)
{
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(), Z, A);
  return fInelasticXsc;
}

/*
inline G4double
G4ComponentGGNuclNuclXsc::GetElementCrossSection(const G4DynamicParticle* dp,
						 G4int Z, const G4Material*)
{
  G4int A = G4lrint(fNist->GetAtomicMassAmu(Z));
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(), Z, A);
  return fInelasticXsc;
}
*/
inline G4double
G4ComponentGGNuclNuclXsc::GetZandACrossSection(const G4DynamicParticle* dp,
					       G4int Z, G4int A)
{
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(), Z, A);
  return fInelasticXsc;
}

inline G4double 
G4ComponentGGNuclNuclXsc::GetCoulombBarier(const G4DynamicParticle* dp, 
					   G4double Z, G4double A, 
					   G4double pR, G4double tR)
{
  return ComputeCoulombBarier(dp->GetDefinition(), dp->GetKineticEnergy(),
                              G4lrint(Z), G4lrint(A), pR, tR);
}

#endif
