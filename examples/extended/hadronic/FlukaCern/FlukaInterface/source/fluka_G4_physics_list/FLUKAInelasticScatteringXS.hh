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
// FLUKA hadron nucleus inelastic XS. 
//
// Author: G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA
#ifndef FLUKA_INELASTIC_SCATTERING_XS_HH
#define FLUKA_INELASTIC_SCATTERING_XS_HH


// G4
#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ElementData.hh"
#include "G4PhysicsVector.hh"


class G4DynamicParticle;
class G4ParticleDefinition;
class G4Element;
class G4Meterial;
class G4Isotope;


class FLUKAInelasticScatteringXS final : public G4VCrossSectionDataSet {

 public: 
  explicit FLUKAInelasticScatteringXS();

  void CrossSectionDescription(std::ostream&) const;
    
  virtual G4bool IsElementApplicable(const G4DynamicParticle*, 
                                     G4int Z,
                                     const G4Material*) override;

  virtual G4bool IsIsoApplicable(const G4DynamicParticle*, 
                                 G4int Z, G4int A,
                                 const G4Element*, const G4Material*) override;

  virtual G4double GetElementCrossSection(const G4DynamicParticle*, 
                                          G4int Z, 
                                          const G4Material*) override;

  virtual G4double GetIsoCrossSection(const G4DynamicParticle*, 
                                      G4int Z, G4int A,
                                      const G4Isotope* iso,
                                      const G4Element* elm,
                                      const G4Material* mat) override;	

  FLUKAInelasticScatteringXS & operator=(const FLUKAInelasticScatteringXS &right) = delete;
  FLUKAInelasticScatteringXS(const FLUKAInelasticScatteringXS&) = delete;
};


#endif
#endif // G4_USE_FLUKA
