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
// File name:     G4VhElectronicStoppingPower
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
// Low energy hadrons/ions electronic stopping power parametrisation
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#ifndef G4VhElectronicStoppingPower_h
#define G4VhElectronicStoppingPower_h 1

#include "G4ios.hh"
#include "globals.hh"

class G4Material;

class G4VhElectronicStoppingPower 
{

public:

  G4VhElectronicStoppingPower();

  virtual ~G4VhElectronicStoppingPower();

  virtual G4double StoppingPower(const G4Material* material,
                                       G4double kineticEnergy) = 0 ;

  virtual G4bool HasMaterial(const G4Material* material) = 0 ;

  virtual G4double ElectronicStoppingPower(G4double z,
                                           G4double kineticEnergy) const = 0;
 
protected:
 
  G4double HeEffChargeSquare(const G4double z, 
                             const G4double kineticEnergyHe) const;
  // This method returns He effective charge square parametrised according to
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985

  G4double GetHeMassAMU() const {return theHeMassAMU;};

private:

  // hide  assignment operator 

    G4VhElectronicStoppingPower(G4VhElectronicStoppingPower &);
    G4VhElectronicStoppingPower & operator =
                    (const G4VhElectronicStoppingPower &right);

  const G4double theHeMassAMU;
};

#endif
