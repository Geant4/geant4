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
// 22.01.07 V.Ivanchenko - add cross section interfaces with Z and A
// 05.03.07 V.Ivanchenko - add IfZAApplicable
//


#ifndef G4UPiNuclearCrossSection_h
#define G4UPiNuclearCrossSection_h

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "globals.hh"
#include "G4Threading.hh"

class G4PhysicsTable;

class G4UPiNuclearCrossSection : public G4VCrossSectionDataSet
{
public:
  
  explicit G4UPiNuclearCrossSection();

  ~G4UPiNuclearCrossSection() override;

  G4bool IsElementApplicable(const G4DynamicParticle* aParticle, 
			     G4int Z, const G4Material*) final;

  G4double GetElasticCrossSection(const G4DynamicParticle* aParticle, 
				  G4int Z, G4int A) const;

  G4double GetInelasticCrossSection(const G4DynamicParticle* aParticle, 
				    G4int Z, G4int A) const;

  void BuildPhysicsTable(const G4ParticleDefinition&) final;

  void DumpPhysicsTable(const G4ParticleDefinition&) final;

  void CrossSectionDescription(std::ostream&) const final;
  
private:

  G4double Interpolate(G4int Z, G4int A, G4double ekin, G4PhysicsTable*) const;

  void AddDataSet(const G4String& p, const G4double* tot, 
		  const G4double* in, const G4double* e, G4int n); 

  void LoadData();

  const G4ParticleDefinition* piPlus;
  const G4ParticleDefinition* piMinus;

  static const G4int NZ = 16;
  static G4int theZ[NZ];

  static G4double theA[NZ];
  static G4double APower[93];

  static G4PhysicsTable* piPlusElastic;
  static G4PhysicsTable* piPlusInelastic;
  static G4PhysicsTable* piMinusElastic;
  static G4PhysicsTable* piMinusInelastic;

  G4double aPower;
  G4double elow;
  G4double elowest;

  G4bool isMaster;

#ifdef G4MULTITHREADED
  static G4Mutex pionUXSMutex;
#endif
};

#endif
