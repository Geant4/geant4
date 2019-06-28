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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BGGNucleonElasticXS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.03.2007
// Modifications:
//
//
// Class Description:
//
// Wrapper of proton and neutron elastic cross-sections using Barashenkov 
// parametersation below 100 GeV and Glauber-Gribov model above
//
// -------------------------------------------------------------------
//

#ifndef G4BGGNucleonElasticXS_h
#define G4BGGNucleonElasticXS_h 1

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"
#include "G4Threading.hh"

class G4ComponentGGHadronNucleusXsc;
class G4NucleonNuclearCrossSection;
class G4HadronNucleonXsc;
class G4Material;
class G4Element;
class G4Isotope;

class G4BGGNucleonElasticXS : public G4VCrossSectionDataSet
{
public:

  explicit G4BGGNucleonElasticXS (const G4ParticleDefinition*);

  ~G4BGGNucleonElasticXS() final;
   
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z, 
			     const G4Material* mat) final;

  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,  
			 const G4Element* elm,
			 const G4Material* mat) final;

  G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z,
				  const G4Material* mat) final;

  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,  
			      const G4Isotope* iso=nullptr,
			      const G4Element* elm=nullptr,
			      const G4Material* mat=nullptr) final;

  void BuildPhysicsTable(const G4ParticleDefinition&) final;

  void CrossSectionDescription(std::ostream&) const final;

private:

  G4double CoulombFactor(G4double kinEnergy, G4int Z);

  G4BGGNucleonElasticXS & operator=(const G4BGGNucleonElasticXS &right);
  G4BGGNucleonElasticXS(const G4BGGNucleonElasticXS&);

  G4double fGlauberEnergy;  
  G4double fLowEnergy;  

  static G4double theGlauberFac[93];
  static G4double theCoulombFac[93];
  static G4int    theA[93];

  const G4ParticleDefinition*     particle;
  const G4ParticleDefinition*     theProton;
  G4ComponentGGHadronNucleusXsc*  fGlauber;
  G4NucleonNuclearCrossSection*   fNucleon;
  G4HadronNucleonXsc*             fHadron;
  G4bool                          isProton;
  G4bool                          isMaster;

#ifdef G4MULTITHREADED
  static G4Mutex nucleonElasticXSMutex;
#endif
};

#endif
