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
// File name:     G4hICRU49He
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
// Electronic stopping power parameterised according to
// ICRU Report N49, 1993. J.F. Ziegler model for He ion.
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4hICRU49He_h
#define G4hICRU49He_h 1

#include "G4VhElectronicStoppingPower.hh"

class G4hICRU49He : public G4VhElectronicStoppingPower
{

public:

  G4hICRU49He();

  ~G4hICRU49He();

  G4bool HasMaterial(const G4Material* material);

  G4double StoppingPower(const G4Material* material,
                               G4double kineticEnergy) ;

  G4double ElectronicStoppingPower(G4double z,
                                   G4double kineticEnergy) const;
 
protected:

private:

  void SetMoleculaNumber(G4int number) {iMolecula = number;};

  G4double rateMass;         // HeMassAMU/ProtonMassAMU
  G4int iMolecula;           // index in the molecula's table
 
};

#endif
