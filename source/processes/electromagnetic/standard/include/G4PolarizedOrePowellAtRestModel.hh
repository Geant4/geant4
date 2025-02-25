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
// GEANT4 Class header file
//
// File name:  G4PolarizedOrePowellAtRestModel
//
// Author:     I.Semeniouk & D.Bernard
//
// Creation date: 26 July 2024
//
// Class Description:
//
// Polarized Ore & Powell AtRest positron 3-gamma annihilation model.
// Electron of media is assumed to have zero kinetic energy.
// Ortho - Positronium three-photon annihilation by Ore and
// Powell, Phys. Rev. 75 (1949).
// Polarization based on R. M. Drisko, Phys. Rev. 102 (1542).
//
// -------------------------------------------------------------------
//

#ifndef G4PolarizedOrePowellAtRestModel_h
#define G4PolarizedOrePowellAtRestModel_h 1

#include "G4VPositronAtRestModel.hh"

class G4PolarizedOrePowellAtRestModel : public G4VPositronAtRestModel
{

public:

  G4PolarizedOrePowellAtRestModel();

  ~G4PolarizedOrePowellAtRestModel() override = default;

  void SampleSecondaries(std::vector<G4DynamicParticle*>& secParticles,
			 G4double&, const G4Material*) const override;

  void PrintGeneratorInformation() const override;

  G4PolarizedOrePowellAtRestModel& operator=
  (const  G4PolarizedOrePowellAtRestModel& right) = delete;
  G4PolarizedOrePowellAtRestModel(const G4PolarizedOrePowellAtRestModel&) = delete;
};

#endif
