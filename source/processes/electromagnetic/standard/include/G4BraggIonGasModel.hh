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
// GEANT4 Class header file
//
//
// File name:     G4BraggIonGasModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 20.05.2010
//
// Modifications:
//
// Class Description:
//
// Implementation of energy loss using and delta-electron production
// by heavy positive partially stripted ion using G4BraggModel. Effective 
// charge of ion is not used, dynamic charge should be sampled by a 
// dedicated charge exchange process
//

// -------------------------------------------------------------------
//

#ifndef G4BraggIonGasModel_h
#define G4BraggIonGasModel_h 1

#include "G4BraggModel.hh"

class G4BraggIonGasModel : public G4BraggModel
{

public:

  explicit G4BraggIonGasModel(const G4ParticleDefinition* p = nullptr,
			      const G4String& nam = "BraggIonGas");

  ~G4BraggIonGasModel() override;

  // Access ion effective charge square ratio to unit charge
  G4double ChargeSquareRatio(const G4Track&) final;

  // Access ion effective charge 
  G4double GetParticleCharge(const G4ParticleDefinition*,
			     const G4Material* mat,
			     G4double kineticEnergy) final;

  // hide assignment operator
  G4BraggIonGasModel & operator=(const  G4BraggIonGasModel &right) = delete;
  G4BraggIonGasModel(const  G4BraggIonGasModel&) = delete;

private:

  G4double currentCharge;
};

#endif
