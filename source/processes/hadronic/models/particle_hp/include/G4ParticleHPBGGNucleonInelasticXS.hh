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
// File name:     G4ParticleHPBGGNucleonInelasticXS
//
// Author:        KOI, Tatsumi
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

// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPBGGNucleonInelasticXS_h
#define G4ParticleHPBGGNucleonInelasticXS_h

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"
#include "G4BGGNucleonInelasticXS.hh"
/*
class G4GlauberGribovCrossSection;
class G4NucleonNuclearCrossSection;
class G4HadronNucleonXsc;
class G4HadronInelasticDataSet;
*/
class G4Material;
class G4Element;
class G4Isotope;

class G4ParticleHPBGGNucleonInelasticXS : public G4BGGNucleonInelasticXS 
{
public:

  G4ParticleHPBGGNucleonInelasticXS (const G4ParticleDefinition*);

  virtual ~G4ParticleHPBGGNucleonInelasticXS();
   
  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z,
			     const G4Material* mat = 0);

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,  
			 const G4Element* elm = 0,
			 const G4Material* mat = 0);

private:
  G4double fLowEnergyLimitForHPN; 
};

#endif
