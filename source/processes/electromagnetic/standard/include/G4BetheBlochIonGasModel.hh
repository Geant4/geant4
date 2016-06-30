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
// $Id: G4BetheBlochIonGasModel.hh 96934 2016-05-18 09:10:41Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BetheBlochIonGasModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 20.05.2010
//
// Modifications:
//
// Class Description:
//
// Implementation of Bethe-Bloch model energy loss and delta-electron production
// by heavy positive partially stripted ion. Effective charge is not used,
// dynamic charge should be sampled by a dedicated charge exchange process

// -------------------------------------------------------------------
//

#ifndef G4BetheBlochIonGasModel_h
#define G4BetheBlochIonGasModel_h 1

#include "G4BetheBlochModel.hh"

class G4BetheBlochIonGasModel : public G4BetheBlochModel
{

public:

  explicit G4BetheBlochIonGasModel(const G4ParticleDefinition* p = 0,
				   const G4String& nam = "BetheBlochGasIon");

  virtual ~G4BetheBlochIonGasModel();

  virtual G4double ChargeSquareRatio(const G4Track& track) final;

  virtual G4double GetParticleCharge(const G4ParticleDefinition* p,
				     const G4Material* mat,
				     G4double kineticEnergy) final;

private:

  // hide assignment operator
  G4BetheBlochIonGasModel & operator=
  (const  G4BetheBlochIonGasModel &right) = delete;
  G4BetheBlochIonGasModel(const  G4BetheBlochIonGasModel&) = delete;

  G4double currentCharge;
};

#endif
