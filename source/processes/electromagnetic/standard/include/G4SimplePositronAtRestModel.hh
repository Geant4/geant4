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
// File name:  G4SimplePositronAtRestModel
//
// Author:     Vladimir Ivanchenko
// 
// Creation date: 14 May 2024
//
// Class Description: 
//
// Simple AtRest positron 2-gamma annihilation model.
// Electron of media is assumed to have zero kinetic energy.
// Polarisation of gamma according to M.H.L.Pryce and J.C.Ward, 
// Nature 4065 (1947) 435.
//
// -------------------------------------------------------------------
//

#ifndef G4SimplePositronAtRestModel_h
#define G4SimplePositronAtRestModel_h 1

#include "G4VPositronAtRestModel.hh"

class G4SimplePositronAtRestModel : public G4VPositronAtRestModel
{

public:

  G4SimplePositronAtRestModel();

  ~G4SimplePositronAtRestModel() override = default;

  void SampleSecondaries(std::vector<G4DynamicParticle*>& secParticles,
			 G4double&, const G4Material*) const override;
  
  void PrintGeneratorInformation() const override;

  G4SimplePositronAtRestModel& operator=
  (const  G4SimplePositronAtRestModel& right) = delete;
  G4SimplePositronAtRestModel(const G4SimplePositronAtRestModel&) = delete;
};

#endif
