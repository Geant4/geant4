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
// Virtual class to provide the interface between G4hLowEnergyIonisation
// and a model of energy loss of low energy proton or alpha particle.
// Further documentation available from http://www.ge.infn.it/geant4/lowE

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
