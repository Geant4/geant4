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
#ifndef G4DiscreteScatteringModel_HH
#define G4DiscreteScatteringModel_HH

#include "G4VEmModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementData.hh"
#include "G4PhysicsVector.hh"
#include "G4ParticleChangeForGamma.hh"

#include <vector>

typedef std::vector<G4double> G4PVDataVector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DiscreteScatteringModel : public G4VEmModel
{

public:
  
  G4DiscreteScatteringModel(G4int iNumAngles=1);

  virtual ~G4DiscreteScatteringModel();
  
  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);
  
  void ReadData(G4int, const G4String & argFileName);
          
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                G4double E, 
                G4double Z, 
                G4double A = 0., 
                G4double cutEnergy = 0.0,
                G4double maxEnergy = DBL_MAX);
                  
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*, 
                                     const G4MaterialCutsCouple*, 
                                     const G4DynamicParticle*, 
                                     G4double tmin=0.0, 
                                     G4double maxEnergy=DBL_MAX);
  
  inline void SetNumberOfAngles(G4int N)        { fNumAngles=N; };
  inline void SetAnalog(const G4String& model)  { fAnalogModel=model; };
      
private:

  G4ThreeVector GetNewDirection(G4double z1);

  static G4ElementData*     fCdf;
  static G4ElementData*     fTcs;
  G4PVDataVector            fGrid;
  G4ParticleChangeForGamma* fParticleChange; 
  G4String                  fAnalogModel;
  G4int                     fNumAngles;
  G4double                  fLowEnergyLimit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


