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
// $Id: G4LivermoreGammaConversionModelRC.hh 94854 2015-12-11 13:54:09Z gcosmo $
//
// Authors: G.Depaola & F.Longo
//

#ifndef G4LivermoreGammaConversionModelRC_h
#define G4LivermoreGammaConversionModelRC_h 1

#include "G4VEmModel.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4ProductionCutsTable.hh"

class G4ParticleChangeForGamma;

class G4LivermoreGammaConversionModelRC : public G4VEmModel
{

public:

  G4LivermoreGammaConversionModelRC(const G4ParticleDefinition* p = 0, 
				    const G4String& nam = "LivermoreGammaConversionRC_1");
  
  virtual ~G4LivermoreGammaConversionModelRC();
  
  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual void InitialiseLocal(const G4ParticleDefinition*, 
			       G4VEmModel* masterModel);
  
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

  virtual G4double MinPrimaryEnergy(const G4Material*,
				    const G4ParticleDefinition*,
				    G4double); 
  
private:
  
  G4ParticleChangeForGamma* fParticleChange;

  G4double lowEnergyLimit;  
  
  G4bool isInitialised;
  G4int verboseLevel;

  void ReadData(size_t Z, const char* path = 0);
  
  G4double ScreenFunction1(G4double screenVariable);
  G4double ScreenFunction2(G4double screenVariable);
  
  G4double smallEnergy;
  //  G4double Psi, Phi;

  
  G4LivermoreGammaConversionModelRC & operator=(const  G4LivermoreGammaConversionModelRC &right);
  G4LivermoreGammaConversionModelRC(const  G4LivermoreGammaConversionModelRC&);

  static G4int maxZ;
  static G4LPhysicsFreeVector* data[100]; // 100 because Z range is 1-99


 G4double fbeta (G4double x);
 G4double Dilog (G4double x);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
