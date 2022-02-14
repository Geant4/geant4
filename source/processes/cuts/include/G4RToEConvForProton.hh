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
// G4RToEConvForProton
//
// Class description:
//
// This class is a Range to Energy Converter for proton.

// Author: H.Kurashige, 05 October 2002 - First implementation
// --------------------------------------------------------------------
#ifndef G4RToEConvForProton_hh
#define G4RToEConvForProton_hh 1

#include "globals.hh"
#include "G4VRangeToEnergyConverter.hh"

class G4RToEConvForProton : public G4VRangeToEnergyConverter
{
public:

  explicit G4RToEConvForProton();

  virtual ~G4RToEConvForProton();

  // Calculate energy cut from given range cut for the material
  G4double Convert(const G4double rangeCut, const G4Material* material) final;

protected:

  G4double ComputeValue(const G4int Z, const G4double kinEnergy) final;

};

#endif
