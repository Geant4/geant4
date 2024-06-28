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
// File name:  G4SimplePsAtRestModel
//
// Author:     I.Semeniouk & D.Bernard
//
// Creation date: 4 Juin 2024
//
// Class Description:
//
// Simple switcher between positron 2-gamma annihilation at rest model and
// positron 2-gamma annihilation at rest model.
//
// -------------------------------------------------------------------
//

#ifndef G4SimplePsAtRestModel_h
#define G4SimplePsAtRestModel_h 1

#include "G4VPositronAtRestModel.hh"

class G4SimplePositronAtRestModel;
class G4OrePowellAtRestModel;

class G4SimplePsAtRestModel : public G4VPositronAtRestModel
{

public:

  G4SimplePsAtRestModel();

  ~G4SimplePsAtRestModel() override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>& secParticles,
			 G4double&, const G4Material*) const override;

  void PrintGeneratorInformation() const override;

  G4SimplePsAtRestModel& operator=
  (const  G4SimplePsAtRestModel& right) = delete;
  G4SimplePsAtRestModel(const G4SimplePsAtRestModel&) = delete;

private:

  G4double f3gFranction;
  G4SimplePositronAtRestModel *model2g;
  G4OrePowellAtRestModel  *model3g;
};

#endif
