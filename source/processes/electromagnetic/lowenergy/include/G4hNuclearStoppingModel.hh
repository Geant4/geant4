//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4hNuclearStoppingModel
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
// Low energy hadrons/ions ionisation parameterisation
// Further documentation available from http://www.ge.infn.it/geant4/lowE/

// -------------------------------------------------------------------
//

#ifndef G4hNuclearStoppingModel_h
#define G4hNuclearStoppingModel_h 1

#include "G4VLowEnergyModel.hh"
#include "G4VhNuclearStoppingPower.hh"

class G4hNuclearStoppingModel : public G4VLowEnergyModel
{

public:

  G4hNuclearStoppingModel(const G4String& name);

  ~G4hNuclearStoppingModel() ;

  G4double TheValue(const G4DynamicParticle* particle,
	            const G4Material* material);

  G4double TheValue(const G4ParticleDefinition* aParticle,
	            const G4Material* material,
                          G4double kineticEnergy);

  G4double HighEnergyLimit(const G4ParticleDefinition* aParticle,
                           const G4Material* material) const
                          {return highEnergyLimit;};
 
  G4double LowEnergyLimit(const G4ParticleDefinition* aParticle,
                          const G4Material* material) const
                          {return lowEnergyLimit;};

  G4double HighEnergyLimit(const G4ParticleDefinition* aParticle) const
                          {return highEnergyLimit;};
 
  G4double LowEnergyLimit(const G4ParticleDefinition* aParticle) const
                          {return lowEnergyLimit;};
 
  G4bool IsInCharge(const G4DynamicParticle* particle,
		    const G4Material* material) const
                          {return true;};

  G4bool IsInCharge(const G4ParticleDefinition* aParticle,
		    const G4Material* material) const
                          {return true;};

  void SetNuclearStoppingFluctuationsOn() 
               {nStopingPowerTable->SetNuclearStoppingFluctuationsOn();}; 

  void SetNuclearStoppingFluctuationsOff() 
               {nStopingPowerTable->SetNuclearStoppingFluctuationsOff();}; 


protected:

private:

  // hide  assignment operator 
  G4hNuclearStoppingModel(G4hNuclearStoppingModel &);
  G4hNuclearStoppingModel & operator=(const G4hNuclearStoppingModel &right);

  void InitializeMe();

  G4double StoppingPower(const G4Material* material,
                               G4double kineticEnergy,
                               G4double z1,
                               G4double m1) const;

  // Pointer to the parametrisation class
  G4VhNuclearStoppingPower* nStopingPowerTable;

  G4double factorPDG2AMU;    // Factor to convert PDG mass unit
                             // into AMU

  G4double theZieglerFactor; // Factor to convert the Stopping Power 
                             // unit [ev/(10^15 atoms/cm^2]
                             // into the Geant4 dE/dx unit
  G4String modelName;
  G4double lowEnergyLimit;
  G4double highEnergyLimit;
};

#endif
