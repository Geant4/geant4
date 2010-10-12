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
// Calculation of hadron elastic cross-sections using selected set of 
// original componetns
//
// 12.08.06 V.Ivanchenko - first implementation
// 22.01.07 V.Ivanchenko - add GetIsoZACrossSection
// 05.03.07 V.Ivanchenko - use G4NucleonNuclearCrossSection
// 06.03.07 V.Ivanchenko - add Initialise function
//
//

#ifndef G4UElasticCrossSection_h
#define G4UElasticCrossSection_h

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4GlauberGribovCrossSection;
class G4NucleonNuclearCrossSection;
class G4UPiNuclearCrossSection;
class G4HadronCrossSections;
class G4ParticleDefinition;

class G4UElasticCrossSection : public G4VCrossSectionDataSet
{
public:

  G4UElasticCrossSection (const G4ParticleDefinition* p=0);

  virtual ~G4UElasticCrossSection ();
   
  virtual
  G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A);

  virtual
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, G4double aTemperature = 0.);

  virtual
  G4double GetZandACrossSection(const G4DynamicParticle*, G4int /*Z*/,
                                G4int /*A*/, G4double aTemperature = 0.);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&); 

private:

  void Initialise(const G4ParticleDefinition*);

  G4bool   hasGlauber;
  G4double thEnergy;  // threshold of Glauber model
  G4double theFac[93];

  G4GlauberGribovCrossSection*    fGlauber;
  G4NucleonNuclearCrossSection*   fNucleon;
  G4UPiNuclearCrossSection*       fUPi;
  G4HadronCrossSections*          fGheisha;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

};

#endif
