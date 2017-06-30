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
// Author: G.Depaola & F.Longo

#ifndef G4BoldyshevTripletModel_h
#define G4BoldyshevTripletModel_h 1

#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Log.hh"

class G4BoldyshevTripletModel : public G4VEmModel
{

public:
  
  explicit G4BoldyshevTripletModel(const G4ParticleDefinition* p = nullptr, 
				   const G4String& nam = "BoldyshevTripletConversion");
  
  virtual ~G4BoldyshevTripletModel();
  
  virtual void Initialise(const G4ParticleDefinition*, 
                          const G4DataVector&);

  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A=0., 
                                      G4double cut=0.,
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

  void ReadData(size_t Z, const char* path = nullptr);

  G4double ScreenFunction1(G4double screenVariable);
  G4double ScreenFunction2(G4double screenVariable);

  G4BoldyshevTripletModel & operator=(const  G4BoldyshevTripletModel &right) = delete;
  G4BoldyshevTripletModel(const  G4BoldyshevTripletModel&) = delete;

  G4int verboseLevel;

  G4double lowEnergyLimit;  
  G4double smallEnergy;
  G4double energyThreshold;
  G4double momentumThreshold_c;
  G4double xb;
  G4double xn;
  
  static G4int maxZ;
  static G4LPhysicsFreeVector* data[100]; // 100 because Z range is 1-99
  
  G4ParticleChangeForGamma* fParticleChange;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
