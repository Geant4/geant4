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
// $Id: G4BGGNucleonInelasticXS.hh 93682 2015-10-28 10:09:49Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BGGNucleonInelasticXS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.03.2007
// Modifications:
//
//
// Class Description:
//
// Wrapper of proton and neutron inelastic cross-sections using Barashenkov 
// parametersation below 100 GeV and Glauber-Gribov model above
//
// -------------------------------------------------------------------
//

#ifndef G4BGGNucleonInelasticXS_h
#define G4BGGNucleonInelasticXS_h

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"

class G4ComponentGGHadronNucleusXsc;
class G4NucleonNuclearCrossSection;
class G4HadronNucleonXsc;
class G4ComponentSAIDTotalXS;
class G4Material;
class G4Element;
class G4Isotope;

class G4BGGNucleonInelasticXS : public G4VCrossSectionDataSet
{
public:

  G4BGGNucleonInelasticXS (const G4ParticleDefinition*);

  virtual ~G4BGGNucleonInelasticXS();
   
  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z,
			     const G4Material* mat = 0);

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,  
			 const G4Element* elm = 0,
			 const G4Material* mat = 0);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z,
				  const G4Material* mat = 0);

  virtual
  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,  
			      const G4Isotope* iso = 0,
			      const G4Element* elm = 0,
			      const G4Material* mat = 0);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual void CrossSectionDescription(std::ostream&) const;

  inline void SetLowestCrossSection(G4double val);

private:

  G4double CoulombFactor(G4double kinEnergy, G4int Z);

  G4BGGNucleonInelasticXS & operator=(const G4BGGNucleonInelasticXS &right);
  G4BGGNucleonInelasticXS(const G4BGGNucleonInelasticXS&);

  G4double fGlauberEnergy;  
  G4double fLowEnergy;  
  G4double fHighEnergy;  
  G4double fSAIDHighEnergyLimit;
  G4double fLowestXSection;
  G4double theGlauberFac[93];
  G4double theCoulombFac[93];
  G4int    theA[93];

  const G4ParticleDefinition*     particle;
  const G4ParticleDefinition*     theProton;

  G4ComponentGGHadronNucleusXsc*  fGlauber;
  G4NucleonNuclearCrossSection*   fNucleon;
  G4HadronNucleonXsc*             fHadron;
  G4ComponentSAIDTotalXS*         fSAID;
  G4bool                          isProton;
  G4bool                          isInitialized;
};

inline void G4BGGNucleonInelasticXS::SetLowestCrossSection(G4double val)
{
  fLowestXSection = val;
}

#endif
