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
// GEANT4 Class file
//
//
// File name:     G4VLowEnergyModel
//
// Author:        Maria Grazia Pia (MariaGrazia.Pia@ge.infn.it)
// 
// Creation date: 7 May 2000
//
// Modifications: 
// 22/05/2000  MGP          Version compliant with design
// 20/07/2000  V.Ivanchenko First implementation
//
// Class Description: 
//
// Abstract class for Low Energy Electromagnetic models
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// -------------------------------------------------------------------
//

#ifndef G4VLowEnergyModel_h
#define G4VLowEnergyModel_h 1

#include "G4ios.hh"
#include "globals.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4Material;

class G4VLowEnergyModel 
{
public:
  explicit G4VLowEnergyModel(const G4String& name);
  virtual ~G4VLowEnergyModel();

  virtual G4double TheValue(const G4DynamicParticle* particle,
			    const G4Material* material)  = 0;

  virtual G4double TheValue(const G4ParticleDefinition* aParticle,
			    const G4Material* material,
                                  G4double kineticEnergy) = 0;

  virtual G4double HighEnergyLimit(const G4ParticleDefinition* aParticle,
                                   const G4Material* material) const = 0;
 
  virtual G4double LowEnergyLimit(const G4ParticleDefinition* aParticle,
                                  const G4Material* material) const = 0;

  virtual G4double HighEnergyLimit(const G4ParticleDefinition* aParticle)
                                  const = 0;
 
  virtual G4double LowEnergyLimit(const G4ParticleDefinition* aParticle) 
                                  const = 0;
 
  virtual G4bool IsInCharge(const G4DynamicParticle* particle,
			    const G4Material* material) const = 0;
 
  virtual G4bool IsInCharge(const G4ParticleDefinition* aParticle,
			    const G4Material* material) const = 0;

  G4VLowEnergyModel & operator=(const  G4VLowEnergyModel &right) = delete;
  G4VLowEnergyModel(const  G4VLowEnergyModel&) = delete;
};
#endif

