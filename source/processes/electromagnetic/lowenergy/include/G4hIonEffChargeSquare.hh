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
// File name:     G4hIonEffChargeSquare
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 20 July 2000
//
// Modifications: 
// 20/07/2000  V.Ivanchenko First implementation
//
// Class Description: 
//
// Ion effective charge model
// J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
// Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4hIonEffChargeSquare_h
#define G4hIonEffChargeSquare_h 1

#include "globals.hh"
#include "G4VLowEnergyModel.hh"

class G4Material;
class G4ParticleDefinition;
class G4DynamicParticle;

class G4hIonEffChargeSquare : public G4VLowEnergyModel
{

public:

  G4hIonEffChargeSquare(const G4String& name) ;

  ~G4hIonEffChargeSquare() ;

  G4double TheValue(const G4DynamicParticle* particle,
	       	          const G4Material* material);

  G4double TheValue(const G4ParticleDefinition* aParticle,
       		          const G4Material* material,
                                G4double kineticEnergy);

  G4double HighEnergyLimit(const G4ParticleDefinition* aParticle,
                           const G4Material* material) const;

  G4double LowEnergyLimit(const G4ParticleDefinition* aParticle,
                          const G4Material* material) const;
 
  G4double HighEnergyLimit(const G4ParticleDefinition* aParticle) const;

  G4double LowEnergyLimit(const G4ParticleDefinition* aParticle) const;
 
  G4bool IsInCharge(const G4DynamicParticle* particle,
		    const G4Material* material) const;

  G4bool IsInCharge(const G4ParticleDefinition* aParticle,
		    const G4Material* material) const;

protected:

private:

  G4double IonEffChargeSquare(const G4Material* material, 
                                    G4double kineticEnergy,
                                    G4double particleMass,
                                    G4double ionCharge) const;
  // This method returns ion effective charge square parametrised 

  const G4double theHeMassAMU;

};

#endif
