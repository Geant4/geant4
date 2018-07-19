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
// $Id: G4BraggNoDeltaModel.hh 97391 2016-06-02 10:08:45Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BraggNoDeltaModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 18.05.2005
//
// Modifications:
//
//
// Class Description:
//
// Implementation of Bethe-Bloch model of energy loss without delta-ray

// -------------------------------------------------------------------
//

#ifndef G4BraggNoDeltaModel_h
#define G4BraggNoDeltaModel_h 1

#include "G4BraggIonModel.hh"

class G4BraggNoDeltaModel : public G4BraggIonModel
{

public:

  explicit G4BraggNoDeltaModel(const G4ParticleDefinition* p = nullptr,
                      const G4String& nam = "BraggNoD");

  virtual ~G4BraggNoDeltaModel();

  virtual G4double ComputeDEDXPerVolume(const G4Material*,
					const G4ParticleDefinition*,
					G4double kineticEnergy,
					G4double cutEnergy) override;

  virtual G4double CrossSectionPerVolume(const G4Material*,
					 const G4ParticleDefinition*,
					 G4double kineticEnergy,
					 G4double cutEnergy,
					 G4double maxEnergy) override;
private:

  // hide assignment operator
  G4BraggNoDeltaModel & operator=(const  G4BraggNoDeltaModel &right) = delete;
  G4BraggNoDeltaModel(const  G4BraggNoDeltaModel&) = delete;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
