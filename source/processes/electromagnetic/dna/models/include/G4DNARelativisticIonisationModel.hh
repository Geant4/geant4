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
// $Id: G4DNARelativisticIonisationModel.hh 90057 2015-05-11 22:25:50Z matkara $
//

#ifndef G4DNARelativisticIonisationModel_h
#define G4DNARelativisticIonisationModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4PhysicsFreeVector.hh"

#include "G4LogLogInterpolation.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4NistManager.hh"

#include "G4DNACrossSectionDataSet.hh"
//#include "G4DNAWaterExcitationStructure.hh"

class G4DNARelativisticIonisationModel: public G4VEmModel
{
public:
  G4DNARelativisticIonisationModel(const G4ParticleDefinition* p = 0,
                           const G4String& nam = "DNARelativisticIonisationModel");

  virtual ~G4DNARelativisticIonisationModel();

  virtual void Initialise(const G4ParticleDefinition*,
                          const G4DataVector& = *(new G4DataVector()));

  virtual G4double CrossSectionPerVolume(const G4Material* material,
                                         const G4ParticleDefinition* p,
                                         G4double ekin,
                                         G4double emin,
                                         G4double emax);

  virtual G4double GetTotalCrossSection  (const G4Material* material,
                                          const G4ParticleDefinition*,
                                          G4double kineticEnergy);
  virtual G4double GetPartialCrossSection(const G4Material* material,
                                          G4int level,
                                          const G4ParticleDefinition*,
                                          G4double kineticEnergy);
  virtual G4double GetDifferentialCrossSection(const G4Material* material,
                                          const G4ParticleDefinition* particle,
                                          G4double kineticEnergy,
                                          G4double secondaryEnergy,
                                          G4int  level);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
                                 G4double maxEnergy);
  virtual void LoadAtomicStates(G4int z, const char *path);
  inline  void SelectStationary       (G4bool input){statCode = input;};
  inline  void SelectFasterComputation(G4bool input){fasterCode = input;}; 
  

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:

  std::vector <G4int   >  iState    [99];
  std::vector <G4int   >  iShell    [99];
  std::vector <G4int   >  iSubShell [99];
  std::vector <G4double>  Nelectrons[99];
  std::vector <G4double>  Ebinding  [99];
  std::vector <G4double>  Ekinetic  [99];

  std::map <G4int, std::vector<G4double> > eVecEZ;

  typedef std::map <G4int, std::map<G4double, std::vector<G4double> > > 
          DeauxDimensionVecMapZ;
  DeauxDimensionVecMapZ eVecEjeEZ;

  typedef std::map <G4int, std::map<G4int,std::map<G4double, 
          std::vector<G4double>  > > > TriDimensionVecMapZ;
  TriDimensionVecMapZ eProbaShellMapZ;

  typedef std::map <G4int, std::map<G4int, std::map<G4double,
          std::map<G4double, G4double> > > > QuadDimensionMapZ;
  QuadDimensionMapZ eDiffCrossSectionDataZ;
  QuadDimensionMapZ eEjectedEnergyDataZ;

  
  G4DNARelativisticIonisationModel & operator
          =(const  G4DNARelativisticIonisationModel &right);
  G4DNARelativisticIonisationModel(const  G4DNARelativisticIonisationModel&);

  G4double     fLowEnergyLimit=0.;
  G4double     fHighEnergyLimit=0.;

  G4bool       isInitialised=false;
  G4bool       statCode=false;
  G4bool       fasterCode=false;
  G4int        verboseLevel=0;

  const std::vector<G4double>*  fMaterialDensity=nullptr;
  const  G4ParticleDefinition*  fParticleDefinition=nullptr;
  G4VAtomDeexcitation*          fAtomDeexcitation=nullptr;

  G4int RandomSelect(const G4Material* material,
                     const G4ParticleDefinition*,
                     G4double kineticEnergy);

  G4double       GetEjectedElectronEnergy   (
                 const G4Material* material,
                 const G4ParticleDefinition* ,
                 G4double energy,
                 G4int shell          );
  G4ThreeVector  GetEjectedElectronDirection(
                 const G4ParticleDefinition* ,
                 G4double energy,G4double secondaryenergy);

  G4double Interpolate     (G4double e1 ,G4double e2 ,G4double e  ,
                            G4double xs1, G4double xs2);
  G4double QuadInterpolator(G4double e11,G4double e12,G4double e21,G4double e22,
                            G4double x11,G4double x12,G4double x21,G4double x22,
                            G4double t1 ,G4double t2 ,G4double t ,G4double e);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
