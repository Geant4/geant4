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

class G4VhEnergyLossModel 
{

public:

  G4VhEnergyLossModel() {};

  virtual ~G4VhEnergyLossModel() {};

  virtual G4double EnergyLoss() const = 0;

  virtual G4double LowEnergyLimit() const = 0;
 
  virtual G4double HighEnergyLimit() const = 0;

  virtual G4bool IsInCharge(G4double energy, 
			    const G4ParticleDefinition* partDef) const = 0;

protected:


private:

}

#endif
