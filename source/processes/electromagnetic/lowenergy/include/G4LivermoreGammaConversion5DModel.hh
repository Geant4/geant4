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
// Author: Zhuxin Li@CENBG
//         11 March 2020
//         on the base of G4LivermoreGammaConversionModel 
//         derives from G4BetheHeitler5DModel               
// -------------------------------------------------------------------

#ifndef G4LivermoreGammaConversion5DModel_h
#define G4LivermoreGammaConversion5DModel_h 1

#include "G4BetheHeitler5DModel.hh"
#include "G4Log.hh"


class G4ParticleChangeForGamma;
class G4LPhysicsFreeVector;
class G4PhysicsLogVector;

class G4LivermoreGammaConversion5DModel : public G4BetheHeitler5DModel
{

public:

  explicit G4LivermoreGammaConversion5DModel(
                      const G4ParticleDefinition* p = nullptr, 
		      const G4String& nam = "Livermore5DConversion");

  virtual ~G4LivermoreGammaConversion5DModel();

  void Initialise(const G4ParticleDefinition*, 
                  const G4DataVector&) override;
  void InitialiseForElement(const G4ParticleDefinition*, 
                                  G4int Z) override;
	
  G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A=0.0, 
                                      G4double cut=0.0,
                                      G4double emax=DBL_MAX) override;


private:

  void ReadData(size_t Z, const char* path = nullptr);

  G4LivermoreGammaConversion5DModel & operator=
  (const  G4LivermoreGammaConversion5DModel &right) = delete;
  G4LivermoreGammaConversion5DModel(const  G4LivermoreGammaConversion5DModel&) = delete;

  static G4double lowEnergyLimit;  
  static G4int verboseLevel;
  static constexpr G4int maxZ =101;
  static G4LPhysicsFreeVector* data[maxZ]; // 101 because Z range is 1-100
   
  G4ParticleChangeForGamma* fParticleChange;
};


#endif
