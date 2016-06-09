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
// $Id: G4BraggNoDeltaModel.hh,v 1.4 2005/10/29 20:09:35 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

  virtual G4double CrossSectionPerVolume(const G4Material*,
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

inline G4double G4BraggNoDeltaModel::CrossSectionPerVolume(
                            const G4Material*,
			    const G4ParticleDefinition*,
			    G4double, G4double, G4double)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
