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
// based on Barashenkov parametrisations of pion data
//
// 16.08.06 V.Ivanchenko - first implementation
//


#ifndef G4UPiNuclearCrossSection_h
#define G4UPiNuclearCrossSection_h

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4DataVector.hh"
#include "globals.hh"
#include <vector>

class G4PhysicsTable;

class G4UPiNuclearCrossSection : public G4VCrossSectionDataSet
{
public:
  
  G4UPiNuclearCrossSection();

  virtual ~G4UPiNuclearCrossSection();

  G4bool IsApplicable(const G4DynamicParticle* aParticle, 
		      const G4Element* anElement);

  G4double GetCrossSection(const G4DynamicParticle* aParticle, 
			   const G4Element* anElement, G4double);

  G4double GetElasticCrossSection(const G4DynamicParticle* aParticle, 
				  const G4Element* anElement);

  G4double GetInelasticCrossSection(const G4DynamicParticle* aParticle, 
				    const G4Element* anElement);

  void BuildPhysicsTable(const G4ParticleDefinition&);

  void DumpPhysicsTable(const G4ParticleDefinition&);
  
private:

  void Initialise();  

  void AddDataSet(const G4String& p, const G4double* tot, 
		  const G4double* in, const G4double* e, G4int n); 

  G4double Interpolate(G4double Z, G4double A, G4double ekin, 
		       G4PhysicsTable*);

  G4int              NZ;
  std::vector<G4int> theZ;
  G4DataVector       theA;
  G4PhysicsTable*    piPlusElastic;
  G4PhysicsTable*    piPlusInelastic;
  G4PhysicsTable*    piMinusElastic;
  G4PhysicsTable*    piMinusInelastic;

  G4double           aPower;
  G4double           elow;
  G4double           elowest;

  const G4ParticleDefinition* piPlus;
  const G4ParticleDefinition* piMinus;
};

inline G4bool G4UPiNuclearCrossSection::IsApplicable(
		     const G4DynamicParticle* part, 
		     const G4Element* elm)
{
  return ((part->GetDefinition() == piMinus || 
	   part->GetDefinition() == piPlus) &&
	  elm->GetZ() > 1.5);
}

inline G4double G4UPiNuclearCrossSection::GetCrossSection(
		      const G4DynamicParticle* dp, 
		      const G4Element* elm, G4double)
{
  return GetInelasticCrossSection(dp, elm);
}

#endif
