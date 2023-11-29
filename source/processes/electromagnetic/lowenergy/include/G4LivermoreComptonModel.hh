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
//
// Author: Alexander Bagulya
//         11 March 2012
//         on the base of G4LivermoreComptonModel
//
// History:
// --------

#ifndef G4LivermoreComptonModel_h
#define G4LivermoreComptonModel_h 1

#include "G4VEmModel.hh"
#include "G4PhysicsFreeVector.hh"

class G4ParticleChangeForGamma;
class G4VAtomDeexcitation;
class G4ShellData;
class G4DopplerProfile;

class G4LivermoreComptonModel : public G4VEmModel
{
public:

  explicit G4LivermoreComptonModel(const G4ParticleDefinition* p = nullptr, 
		          const G4String& nam = "LivermoreCompton");
  ~G4LivermoreComptonModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  void InitialiseLocal(const G4ParticleDefinition*,
		       G4VEmModel* masterModel) override;
  
  void InitialiseForElement(const G4ParticleDefinition*, G4int Z) override;

  G4double ComputeCrossSectionPerAtom( const G4ParticleDefinition*,
				       G4double kinEnergy, 
				       G4double Z, 
				       G4double A=0, 
				       G4double cut=0,
				       G4double emax=DBL_MAX ) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;
  
  G4LivermoreComptonModel & operator=(const  G4LivermoreComptonModel &right) = delete;
  G4LivermoreComptonModel(const  G4LivermoreComptonModel&) = delete;

private:
  void ReadData(const G4int ZZ);
  const G4String& FindDirectoryPath();
  G4double ComputeScatteringFunction(G4double x, G4int Z);

  G4ParticleChangeForGamma* fParticleChange;
  G4VAtomDeexcitation*      fAtomDeexcitation;

  static G4ShellData*       shellData;
  static G4DopplerProfile*  profileData;


  static G4PhysicsFreeVector* data[101]; // 101 because Z range is 1-100
  static const G4double ScatFuncFitParam[101][16];
  static G4String gDataDirectory;

  G4int verboseLevel;
  G4int maxZ = 100;
  G4bool isInitialised = false;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
