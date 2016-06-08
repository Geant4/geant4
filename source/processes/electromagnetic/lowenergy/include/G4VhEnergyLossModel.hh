// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// For information related to this code contact:
// Geant4 Collaboration
//
// File name:     G4VhEnergyLossModel
//
// Author:        Maria Grazia Pia (MariaGrazia.Pia@ge.infn.it)
// 
// Creation date: 7 May 2000
//
// Modifications: 
// 22/05/2000  MGP  Version compliant with design
//
// Class Description: 
//
// Abstract base class for hadron energy loss model
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#ifndef G4VHENERGYLOSSMODEL_HH
#define G4VHENERGYLOSSMODEL_HH

#include "globals.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4Material;

class G4VhEnergyLossModel 
{

public:

  G4VhEnergyLossModel() {};

  virtual ~G4VhEnergyLossModel() {};

  virtual G4double EnergyLoss(const G4DynamicParticle* particle,
			      const G4Material* material) const = 0;

  virtual G4double LowEnergyLimit() const = 0;
 
  virtual G4double HighEnergyLimit() const = 0;

  virtual G4bool IsInCharge(G4double energy, 
			    const G4ParticleDefinition* partDef,
			    const G4Material* material) const = 0;

protected:


private:

};

#endif
