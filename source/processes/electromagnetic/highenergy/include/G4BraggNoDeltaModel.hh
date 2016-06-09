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
// $Id: G4BraggNoDeltaModel.hh,v 1.6 2006/06/29 19:32:16 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
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

  G4BraggNoDeltaModel(const G4ParticleDefinition* p = 0,
                      const G4String& nam = "BraggNoD");

  virtual ~G4BraggNoDeltaModel();

  virtual G4double ComputeDEDXPerVolume(const G4Material*,
					const G4ParticleDefinition*,
					G4double kineticEnergy,
					G4double cutEnergy);

  virtual G4double ComputeCrossSectionPerElectron(
					 const G4ParticleDefinition*,
					 G4double kineticEnergy,
					 G4double cutEnergy,
					 G4double maxEnergy);

private:

  // hide assignment operator
  G4BraggNoDeltaModel & operator=(const  G4BraggNoDeltaModel &right);
  G4BraggNoDeltaModel(const  G4BraggNoDeltaModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4BraggNoDeltaModel::ComputeDEDXPerVolume(
                            const G4Material* material,
			    const G4ParticleDefinition* pd,
                            G4double kinEnergy, G4double)
{
  G4double dedx = G4BraggIonModel::ComputeDEDXPerVolume(material, pd, kinEnergy, DBL_MAX);
  return dedx;
}

inline G4double G4BraggNoDeltaModel::ComputeCrossSectionPerElectron(
			    const G4ParticleDefinition*,
			    G4double, G4double, G4double)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
