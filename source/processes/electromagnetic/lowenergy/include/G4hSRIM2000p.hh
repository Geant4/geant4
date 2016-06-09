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
// File name:     G4hSRIM2000p
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 7 Jan 2003
//
// Modifications: 
// 07/01/2003  V.Ivanchenko First implementation
//
// Class Description: 
//
// Electronic stopping power parametrised according to
// J. F. Ziegler, SRIM-2000 software package from www.SRIM.org
//

// -------------------------------------------------------------------
//

#ifndef G4hSRIM2000p_h
#define G4hSRIM2000p_h 1

#include "globals.hh"
#include "G4VhElectronicStoppingPower.hh"

class G4Material;

class G4hSRIM2000p : public G4VhElectronicStoppingPower
{

public:

  G4hSRIM2000p();

  ~G4hSRIM2000p();

  G4bool HasMaterial(const G4Material* material);

  G4double StoppingPower(const G4Material* material,
                               G4double kineticEnergy);

  G4double ElectronicStoppingPower(G4double z,
                                   G4double kineticEnergy) const;
 
protected:

private:
  const G4double protonMassAMU;
  static G4double a[92][8];

};

#endif
