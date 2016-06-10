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
// $Id: G4LivermorePolarizedPhotoElectricGDModel.hh 66241 2012-12-13 18:34:42Z gunter $
//
// Authors: G.Depaola & F.Longo
//

#ifndef G4LivermorePolarizedPhotoElectricGDModel_h
#define G4LivermorePolarizedPhotoElectricGDModel_h 1

#include "G4VEmModel.hh"
#include "G4ElementData.hh"

class G4ParticleChangeForGamma;
class G4VAtomDeexcitation; 

#include "G4LPhysicsFreeVector.hh"
#include <vector>

class G4LivermorePolarizedPhotoElectricGDModel : public G4VEmModel
{

public:

  G4LivermorePolarizedPhotoElectricGDModel( 
	     const G4String& nam = "LivermorePolarizedPhotoElectric");

  virtual ~G4LivermorePolarizedPhotoElectricGDModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

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

  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z);
  
  inline void SetLimitNumberOfShells(G4int);

protected:

  G4ParticleChangeForGamma* fParticleChange;

private:
  
  void ReadData(G4int Z, const char* path = 0);

  G4LivermorePolarizedPhotoElectricGDModel & operator=(const  G4LivermorePolarizedPhotoElectricGDModel &right);
  G4LivermorePolarizedPhotoElectricGDModel(const  G4LivermorePolarizedPhotoElectricGDModel&);

  G4ParticleDefinition*   theGamma;
  G4ParticleDefinition*   theElectron;

  G4int    verboseLevel;

  G4int                   maxZ;
  G4int                   nShellLimit;
  G4bool                  fDeexcitationActive;
  G4bool                  isInitialised;

  static G4LPhysicsFreeVector*   fCrossSection[99];
  static G4LPhysicsFreeVector*   fCrossSectionLE[99];
  static std::vector<G4double>*  fParam[99];
  static G4int                   fNShells[99];
  static G4int                   fNShellsUsed[99];
  static G4ElementData*          fShellCrossSection;
  static G4Material*             fWater;
  static G4double                fWaterEnergyLimit;
  
  G4VAtomDeexcitation*    fAtomDeexcitation;
  
  G4double                fCurrSection;
  std::vector<G4double>   fSandiaCof;

  // specific methods for polarization -- FL & GD 
  
  G4ThreeVector GetRandomPolarization(G4ThreeVector& direction0); // Random Polarization
  G4ThreeVector GetPerpendicularPolarization(const G4ThreeVector& direction0, const G4ThreeVector& polarization0) const;
  
  G4ThreeVector SetPerpendicularVector(G4ThreeVector& a); // temporary
  G4ThreeVector SetNewPolarization(G4double epsilon, G4double sinSqrTheta, 
				   G4double phi, G4double cosTheta);
  G4double SetPhi(G4double, G4double, G4double);
  G4double SetCosTheta(G4double);
  
  void SystemOfRefChange(G4ThreeVector& direction0, G4ThreeVector& direction1, 
			 G4ThreeVector& polarization0);

};

inline void G4LivermorePolarizedPhotoElectricGDModel::SetLimitNumberOfShells(G4int n)
{
  nShellLimit = n;
}



#endif
