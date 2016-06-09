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
//
//

#ifndef G4UInelasticCrossSection_h
#define G4UInelasticCrossSection_h

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4GlauberGribovCrossSection;
class G4ProtonInelasticCrossSection;
class G4NeutronInelasticCrossSection;
class G4UPiNuclearCrossSection;
class G4HadronCrossSections;
class G4ParticleDefinition;

class G4UInelasticCrossSection : public G4VCrossSectionDataSet
{
public:

  G4UInelasticCrossSection (const G4ParticleDefinition*);

  virtual ~G4UInelasticCrossSection ();
   
  virtual
  G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

  virtual
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, G4double aTemperature);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&); 

private:

  G4int idx;

  G4GlauberGribovCrossSection*    fGlauber;
  G4ProtonInelasticCrossSection*  fGhadProton;
  G4NeutronInelasticCrossSection* fGhadNeutron;
  G4UPiNuclearCrossSection*       fUPi;
  G4HadronCrossSections*          fGheisha;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

};

#endif
