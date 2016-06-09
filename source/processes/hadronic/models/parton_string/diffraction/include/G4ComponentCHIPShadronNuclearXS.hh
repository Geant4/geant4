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
//  Calculation of the total, elastic and inelastic cross-sections
//  of hadron (proton, neutron, pi+, pi-, K+, K-, anti_proton, anti_neutron
//  interactions with nuclei based on CHIPS model
//
//   Created by V. Uzhinsky, 31.05.2011  


#ifndef G4ComponentCHIPShadronNuclearXS_h
#define G4ComponentCHIPShadronNuclearXS_h

#include <CLHEP/Units/PhysicalConstants.h>  // pi, fermi,..

#include "globals.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Nucleus.hh"

#include "G4VComponentCrossSection.hh"

class G4ParticleDefinition;
class G4VQCrossSection;

class G4ComponentCHIPShadronNuclearXS : public G4VComponentCrossSection
{
public:

  G4ComponentCHIPShadronNuclearXS ();
  virtual ~G4ComponentCHIPShadronNuclearXS ();

  virtual
  G4double GetTotalElementCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy, 
				       G4int Z, G4double N);

  virtual
  G4double GetTotalIsotopeCrossSection(const G4ParticleDefinition* aParticle,
				       G4double kinEnergy,
				       G4int Z, G4int N);
  virtual
  G4double GetInelasticElementCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4double N);
  virtual
  G4double GetInelasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					   G4double kinEnergy, 
					   G4int Z, G4int N);

  virtual
  G4double GetElasticElementCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4double N);

  virtual
  G4double GetElasticIsotopeCrossSection(const G4ParticleDefinition* aParticle,
					 G4double kinEnergy, 
					 G4int Z, G4int N);
 
  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&)
  {}

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&) 
  {}

  private:
  G4ComponentCHIPShadronNuclearXS & operator=(const G4ComponentCHIPShadronNuclearXS &right);
  G4ComponentCHIPShadronNuclearXS(const G4ComponentCHIPShadronNuclearXS&);

  const G4double fUpperLimit;
  const G4double fLowerLimit; 

  G4VQCrossSection* PxsManagerEl;
  G4VQCrossSection* PxsManagerInEl;

  G4VQCrossSection* NxsManagerEl;
  G4VQCrossSection* NxsManagerInEl;

  G4VQCrossSection* PBARxsManagerEl;
  G4VQCrossSection* PBARxsManagerInEl;

  G4VQCrossSection* PIPxsManagerEl; 
  G4VQCrossSection* PIPxsManagerInEl; 

  G4VQCrossSection* PIMxsManagerEl; 
  G4VQCrossSection* PIMxsManagerInEl; 

  G4VQCrossSection* KPxsManagerEl;
  G4VQCrossSection* KPxsManagerInEl;

  G4VQCrossSection* KMxsManagerEl;
  G4VQCrossSection* KMxsManagerInEl;
};

#endif
