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
// $Id: G4LivermorePolarizedComptonModel.hh 93359 2015-10-19 13:42:18Z gcosmo $
//
// Authors: G.Depaola & F.Longo
//

#ifndef G4LivermorePolarizedComptonModel_h
#define G4LivermorePolarizedComptonModel_h 1

#include "G4VEmModel.hh"
#include "G4LPhysicsFreeVector.hh"

class G4ParticleChangeForGamma;
class G4VAtomDeexcitation;
class G4ShellData;
class G4DopplerProfile;

#include "G4VEMDataSet.hh"
#include "G4CompositeEMDataSet.hh"


//#include "G4Electron.hh"
//#include "G4CrossSectionHandler.hh"
//#include "G4LogLogInterpolation.hh"
//#include "G4CompositeEMDataSet.hh"
//#include "G4ForceCondition.hh"
//#include "G4Gamma.hh"

class G4LivermorePolarizedComptonModel : public G4VEmModel
{

public:

  G4LivermorePolarizedComptonModel(const G4ParticleDefinition* p = 0, 
		                   const G4String& nam = "LivermorePolarizedCompton");

  virtual ~G4LivermorePolarizedComptonModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual void InitialiseLocal(const G4ParticleDefinition*, G4VEmModel* masterModel);

  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A=0, 
                                      G4double cut=0,
                                      G4double emax=DBL_MAX);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

private:

  G4ParticleChangeForGamma* fParticleChange;
  G4VAtomDeexcitation* fAtomDeexcitation;

  G4bool isInitialised;
  G4int verboseLevel;

  G4LivermorePolarizedComptonModel & operator=(const  G4LivermorePolarizedComptonModel &right);
  G4LivermorePolarizedComptonModel(const  G4LivermorePolarizedComptonModel&);

  // specific methods for polarization 
  
  G4ThreeVector GetRandomPolarization(G4ThreeVector& direction0); // Random Polarization
  G4ThreeVector GetPerpendicularPolarization(const G4ThreeVector& direction0, const G4ThreeVector& polarization0) const;
  G4ThreeVector SetPerpendicularVector(G4ThreeVector& a); // temporary
  G4ThreeVector SetNewPolarization(G4double epsilon, G4double sinSqrTheta, 
				   G4double phi, G4double cosTheta);
  G4double SetPhi(G4double, G4double);
  void SystemOfRefChange(G4ThreeVector& direction0, G4ThreeVector& direction1, 
			 G4ThreeVector& polarization0, G4ThreeVector& polarization1);
  
  // Doppler Broadening
  
  static G4ShellData* shellData;
  static G4DopplerProfile* profileData;

  // Cross Section Handling 
 
  static G4int maxZ;
  static G4LPhysicsFreeVector* data[100];
  void ReadData(size_t Z, const char* path = 0);

  // Scattering function 

  static G4CompositeEMDataSet* scatterFunctionData;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
