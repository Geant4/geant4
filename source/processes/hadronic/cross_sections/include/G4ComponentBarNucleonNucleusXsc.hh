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
// author: Vladimir.Grichine@cern.ch
//
// Implements data from: Barashenkov V.S., Nucleon-Nucleus Cross Section,
// Preprint JINR P2-89-770, p. 12, Dubna 1989 (scanned version from KEK)
// Based on G4NucleonNuclearCrossSection class
//
// Modifications: 16.08.2018 V.Ivanchenko major revision
//

#ifndef G4ComponentBarNucleonNucleusXsc_h
#define G4ComponentBarNucleonNucleusXsc_h


#include "G4VComponentCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"

#include "globals.hh"
#include "G4PiData.hh"

class G4ComponentBarNucleonNucleusXsc : public G4VComponentCrossSection
{

public:
  
  G4ComponentBarNucleonNucleusXsc();
  virtual ~G4ComponentBarNucleonNucleusXsc();

  virtual
  G4double GetTotalIsotopeCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy,
				       G4int Z, G4int );

  virtual
  G4double GetTotalElementCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy, 
				       G4int Z, G4double );

  virtual
  G4double GetInelasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4int );

  virtual
  G4double GetInelasticElementCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4double );

  virtual
  G4double GetElasticElementCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4double );

  virtual
  G4double GetElasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int );

  void ComputeCrossSections(const G4ParticleDefinition* aParticle,
			    G4double kinEnergy, G4int Z);
 
  G4bool IsElementApplicable(const G4DynamicParticle* aParticle, G4int Z);

  virtual void CrossSectionDescription(std::ostream&) const;

  inline G4double GetElementCrossSection(const G4DynamicParticle* aParticle, G4int Z); 
  inline G4double GetElasticCrossSection(const G4DynamicParticle* aParticle, G4int Z);

  inline G4double GetTotalXsc()     { return fTotalXsc;   };
  inline G4double GetElasticXsc()   { return fElasticXsc; };
  inline G4double GetInelasticXsc() { return fInelasticXsc; };
  
private:

  G4double Interpolate(G4int Z1, G4int Z2, G4int Z, G4double x1, G4double x2);

  std::vector< G4int >     theZ;
  std::vector< G4PiData* > thePipData;
  std::vector< G4PiData* > thePimData;

  // cross sections
  G4double fTotalXsc;
  G4double fInelasticXsc;
  G4double fElasticXsc;

  // particles
  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theNeutron;

  static G4double theA[93];
  static G4double A75[93];

};

inline
G4double G4ComponentBarNucleonNucleusXsc::GetElementCrossSection(
         const G4DynamicParticle* dp, G4int Z)
{
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(),Z);
  return fInelasticXsc;
}

inline
G4double G4ComponentBarNucleonNucleusXsc::GetElasticCrossSection(
         const G4DynamicParticle* dp, G4int Z)
{
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(),Z);
  return fElasticXsc;
}

#endif
