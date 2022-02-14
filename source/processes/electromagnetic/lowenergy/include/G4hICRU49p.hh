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
// File name:     G4hICRU49p
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
// Electronic stopping power parametrised according to
// ICRU Report N49, 1993, for protons.
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#ifndef G4hICRU49p_h
#define G4hICRU49p_h 1

#include "globals.hh"
#include "G4VhElectronicStoppingPower.hh"

class G4Material;

class G4hICRU49p : public G4VhElectronicStoppingPower
{
public:

  explicit G4hICRU49p();

  ~G4hICRU49p();

  G4bool HasMaterial(const G4Material* material) override;

  G4double StoppingPower(const G4Material* material,
                               G4double kineticEnergy) override;

  G4double ElectronicStoppingPower(G4double z,
                                   G4double kineticEnergy) const override;
private:
  void SetMoleculaNumber(G4int number) {iMolecula = number;};

  const G4double protonMassAMU;
  G4int iMolecula;               // index in the molecula's table
};
 
#endif
