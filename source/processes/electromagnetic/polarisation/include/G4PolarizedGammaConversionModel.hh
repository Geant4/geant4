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
// $Id: G4PolarizedGammaConversionModel.hh 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PolarizedGammaConversionModel
//
// Author:        Karim Laihem
//
// Creation date: 19.04.2005
//
// Modifications:
// 21-08-06 Modified to work in g4.8.1 framework (A.Schaelicke)
//
// Class Description:
//
// Implementation of gamma convertion to e+e- in the field of a nucleus 
// including polarization transfer

// -------------------------------------------------------------------
//

#ifndef G4PolarizedGammaConversionModel_h
#define G4PolarizedGammaConversionModel_h 1

#include "G4BetheHeitlerModel.hh"
#include "G4PhysicsTable.hh"

class G4ParticleChangeForGamma;
class G4VPolarizedCrossSection;
class G4PolarizedGammaConversionModel : public G4BetheHeitlerModel
{

public:

  explicit G4PolarizedGammaConversionModel(const G4ParticleDefinition* p = nullptr, 
		      const G4String& nam = "polConv");

  virtual ~G4PolarizedGammaConversionModel();
 
  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;

  inline const G4Element* SelectedAtom();

 protected:
  G4VPolarizedCrossSection*           crossSectionCalculator;

};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4Element* G4PolarizedGammaConversionModel::SelectedAtom()
{
  return GetCurrentElement();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
