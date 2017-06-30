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
// $Id: G4VEmModel.hh 104457 2017-05-31 15:52:37Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VEmModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 23-12-02 V.Ivanchenko change interface before move to cut per region
// 24-01-03 Cut per region (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 25-02-03 Add sample theta and displacement (V.Ivanchenko)
// 23-07-03 Replace G4Material by G4MaterialCutCouple in dE/dx and CrossSection
//          calculation (V.Ivanchenko)
// 01-03-04 L.Urban signature changed in SampleCosineTheta 
// 23-04-04 L.urban signature of SampleCosineTheta changed back 
// 17-11-04 Add method CrossSectionPerAtom (V.Ivanchenko)
// 14-03-05 Reduce number of pure virtual methods and make inline part 
//          separate (V.Ivanchenko)
// 24-03-05 Remove IsInCharge and add G4VParticleChange in the constructor (VI)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 15-04-05 optimize internal interface for msc (V.Ivanchenko)
// 08-05-05 A -> N (V.Ivanchenko)
// 25-07-05 Move constructor and destructor to the body (V.Ivanchenko)
// 02-02-06 ComputeCrossSectionPerAtom: default value A=0. (mma)
// 06-02-06 add method ComputeMeanFreePath() (mma)
// 07-03-06 Optimize msc methods (V.Ivanchenko)
// 29-06-06 Add member currentElement and Get/Set methods (V.Ivanchenko)
// 29-10-07 Added SampleScattering (V.Ivanchenko)
// 15-07-08 Reorder class members and improve comments (VI)
// 21-07-08 Added vector of G4ElementSelector and methods to use it (VI)
// 12-09-08 Added methods GetParticleCharge, GetChargeSquareRatio, 
//          CorrectionsAlongStep, ActivateNuclearStopping (VI)
// 16-02-09 Moved implementations of virtual methods to source (VI)
// 07-04-09 Moved msc methods from G4VEmModel to G4VMscModel (VI)
// 13-10-10 Added G4VEmAngularDistribution (VI)
//
// Class Description:
//
// Abstract interface to energy loss models

// -------------------------------------------------------------------
//

#ifndef G4VEmModel_h
#define G4VEmModel_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4Isotope.hh"
#include "G4DataVector.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4VEmAngularDistribution.hh"
#include "G4EmElementSelector.hh"
#include <CLHEP/Random/RandomEngine.h>
#include <vector>

class G4ElementData;
class G4PhysicsTable;
class G4Region;
class G4VParticleChange;
class G4ParticleChangeForLoss;
class G4ParticleChangeForGamma;
class G4Track;
class G4LossTableManager;

class G4VEmModel
{

public:

  explicit G4VEmModel(const G4String& nam);

  virtual ~G4VEmModel();

  //------------------------------------------------------------------------
  // Virtual methods to be implemented for any concrete model
  //------------------------------------------------------------------------

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&) = 0;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin = 0.0,
                                 G4double tmax = DBL_MAX) = 0;

  //------------------------------------------------------------------------
  // Methods for initialisation of MT; may be overwritten if needed
  //------------------------------------------------------------------------

  // initilisation in local thread
  virtual void InitialiseLocal(const G4ParticleDefinition*,
                               G4VEmModel* masterModel);

  // initilisation of a new material at run time
  virtual void InitialiseForMaterial(const G4ParticleDefinition*,
                                     const G4Material*);

  // initilisation of a new element at run time
  virtual void InitialiseForElement(const G4ParticleDefinition*,
                                    G4int Z);

  //------------------------------------------------------------------------
  // Methods with standard implementation; may be overwritten if needed 
  //------------------------------------------------------------------------

  // main method to compute dEdx
  virtual G4double ComputeDEDXPerVolume(const G4Material*,
                                        const G4ParticleDefinition*,
                                        G4double kineticEnergy,
                                        G4double cutEnergy = DBL_MAX);

  // main method to compute cross section per Volume
  virtual G4double CrossSectionPerVolume(const G4Material*,
                                         const G4ParticleDefinition*,
                                         G4double kineticEnergy,
                                         G4double cutEnergy = 0.0,
                                         G4double maxEnergy = DBL_MAX);

  // method to get partial cross section
  virtual G4double GetPartialCrossSection(const G4Material*,
                                          G4int level,
                                          const G4ParticleDefinition*,
                                          G4double kineticEnergy);

  // main method to compute cross section per atom
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                              G4double kinEnergy,
                                              G4double Z,
                                              G4double A = 0., /* amu */
                                              G4double cutEnergy = 0.0,
                                              G4double maxEnergy = DBL_MAX);

  // main method to compute cross section per atomic shell
  virtual G4double ComputeCrossSectionPerShell(const G4ParticleDefinition*,
                                               G4int Z, G4int shellIdx,
                                               G4double kinEnergy,
                                               G4double cutEnergy = 0.0,
                                               G4double maxEnergy = DBL_MAX);

  // Compute effective ion charge square
  virtual G4double ChargeSquareRatio(const G4Track&);

  // Compute effective ion charge square
  virtual G4double GetChargeSquareRatio(const G4ParticleDefinition*,
                                        const G4Material*,
                                        G4double kineticEnergy);

  // Compute ion charge 
  virtual G4double GetParticleCharge(const G4ParticleDefinition*,
                                     const G4Material*,
                                     G4double kineticEnergy);

  // Initialisation for a new track
  virtual void StartTracking(G4Track*);

  // add correction to energy loss and compute non-ionizing energy loss
  virtual void CorrectionsAlongStep(const G4MaterialCutsCouple*,
                                    const G4DynamicParticle*,
                                    G4double& eloss,
                                    G4double& niel,
                                    G4double length);

  // value which may be tabulated (by default cross section)
  virtual G4double Value(const G4MaterialCutsCouple*,
                         const G4ParticleDefinition*,
                         G4double kineticEnergy);

  // threshold for zero value 
  virtual G4double MinPrimaryEnergy(const G4Material*,
                                    const G4ParticleDefinition*,
                                    G4double cut = 0.0);

  // model can define low-energy limit for the cut
  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
                                const G4MaterialCutsCouple*);

  // initilisation at run time for a given material
  virtual void SetupForMaterial(const G4ParticleDefinition*,
                                const G4Material*,
                                G4double kineticEnergy);

  // add a region for the model
  virtual void DefineForRegion(const G4Region*);

  // for automatic documentation
  virtual void ModelDescription(std::ostream& outFile) const; 

protected:

  // initialisation of the ParticleChange for the model
  G4ParticleChangeForLoss* GetParticleChangeForLoss();

  // initialisation of the ParticleChange for the model
  G4ParticleChangeForGamma* GetParticleChangeForGamma();

  // kinematically allowed max kinetic energy of a secondary
  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                                      G4double kineticEnergy);

public:

  //------------------------------------------------------------------------
  // Generic methods common to all models
  //------------------------------------------------------------------------

  // should be called at initialisation to build element selectors
  void InitialiseElementSelectors(const G4ParticleDefinition*,
                                  const G4DataVector&);

  // should be called at initialisation to access element selectors
  inline std::vector<G4EmElementSelector*>* GetElementSelectors();

  // should be called at initialisation to set element selectors
  inline void SetElementSelectors(std::vector<G4EmElementSelector*>*);

  // dEdx per unit length
  virtual inline G4double ComputeDEDX(const G4MaterialCutsCouple*,
                              const G4ParticleDefinition*,
                              G4double kineticEnergy,
                              G4double cutEnergy = DBL_MAX);

  // cross section per volume
  inline G4double CrossSection(const G4MaterialCutsCouple*,
                               const G4ParticleDefinition*,
                               G4double kineticEnergy,
                               G4double cutEnergy = 0.0,
                               G4double maxEnergy = DBL_MAX);

  // compute mean free path via cross section per volume
  inline G4double ComputeMeanFreePath(const G4ParticleDefinition*,
                                      G4double kineticEnergy,
                                      const G4Material*,
                                      G4double cutEnergy = 0.0,
                                      G4double maxEnergy = DBL_MAX);

  // generic cross section per element
  inline G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                             const G4Element*,
                                             G4double kinEnergy,
                                             G4double cutEnergy = 0.0,
                                             G4double maxEnergy = DBL_MAX);

  // select isotope in order to have precise mass of the nucleus
  inline G4int SelectIsotopeNumber(const G4Element*);

  // atom can be selected effitiantly if element selectors are initialised 
  inline const G4Element* SelectRandomAtom(const G4MaterialCutsCouple*,
                                           const G4ParticleDefinition*,
                                           G4double kineticEnergy,
                                           G4double cutEnergy = 0.0,
                                           G4double maxEnergy = DBL_MAX);

  // to select atom cross section per volume is recomputed for each element 
  const G4Element* SelectRandomAtom(const G4Material*,
                                    const G4ParticleDefinition*,
                                    G4double kineticEnergy,
                                    G4double cutEnergy = 0.0,
                                    G4double maxEnergy = DBL_MAX);

  // to select atom if cross section is proportional number of electrons 
  inline G4int SelectRandomAtomNumber(const G4Material*);

  //------------------------------------------------------------------------
  // Get/Set methods
  //------------------------------------------------------------------------

  void SetParticleChange(G4VParticleChange*, G4VEmFluctuationModel* f=nullptr);

  void SetCrossSectionTable(G4PhysicsTable*, G4bool isLocal);

  inline G4ElementData* GetElementData();

  inline G4PhysicsTable* GetCrossSectionTable();

  inline G4VEmFluctuationModel* GetModelOfFluctuations();

  inline G4VEmAngularDistribution* GetAngularDistribution();

  inline G4VEmModel* GetTripletModel();

  inline void SetTripletModel(G4VEmModel*);

  inline void SetAngularDistribution(G4VEmAngularDistribution*);

  inline G4double HighEnergyLimit() const;

  inline G4double LowEnergyLimit() const;

  inline G4double HighEnergyActivationLimit() const;

  inline G4double LowEnergyActivationLimit() const;

  inline G4double PolarAngleLimit() const;

  inline G4double SecondaryThreshold() const;

  inline G4bool LPMFlag() const;

  inline G4bool DeexcitationFlag() const;

  inline G4bool ForceBuildTableFlag() const;

  inline G4bool UseAngularGeneratorFlag() const;

  inline void SetAngularGeneratorFlag(G4bool);

  inline void SetHighEnergyLimit(G4double);

  inline void SetLowEnergyLimit(G4double);

  inline void SetActivationHighEnergyLimit(G4double);

  inline void SetActivationLowEnergyLimit(G4double);

  inline G4bool IsActive(G4double kinEnergy);

  inline void SetPolarAngleLimit(G4double);

  inline void SetSecondaryThreshold(G4double);

  inline void SetLPMFlag(G4bool val);

  inline void SetDeexcitationFlag(G4bool val);

  inline void SetForceBuildTable(G4bool val);

  inline void SetFluctuationFlag(G4bool val);

  inline void SetMasterThread(G4bool val);

  inline G4bool IsMaster() const;

  inline G4double MaxSecondaryKinEnergy(const G4DynamicParticle* dynParticle);

  inline const G4String& GetName() const;

  inline void SetCurrentCouple(const G4MaterialCutsCouple*);

  inline const G4Element* GetCurrentElement() const;

  inline const G4Isotope* GetCurrentIsotope() const;

  inline G4bool IsLocked() const;

  inline void SetLocked(G4bool);

protected:

  inline const G4MaterialCutsCouple* CurrentCouple() const;

  inline void SetCurrentElement(const G4Element*);

private:

  //  hide assignment operator
  G4VEmModel & operator=(const  G4VEmModel &right) = delete;
  G4VEmModel(const  G4VEmModel&) = delete;

  // ======== Parameters of the class fixed at construction =========
 
  G4VEmFluctuationModel* flucModel;
  G4VEmAngularDistribution* anglModel;
  const G4String   name;

  // ======== Parameters of the class fixed at initialisation =======

  G4double        lowLimit;
  G4double        highLimit;
  G4double        eMinActive;
  G4double        eMaxActive;
  G4double        polarAngleLimit;
  G4double        secondaryThreshold;
  G4bool          theLPMflag;
  G4bool          flagDeexcitation;
  G4bool          flagForceBuildTable;
  G4bool          isMaster;

  G4bool          localTable;
  G4bool          localElmSelectors;
  G4bool          useAngularGenerator;
  G4bool          isLocked;
  G4int           nSelectors;
  std::vector<G4EmElementSelector*>* elmSelectors;
  G4LossTableManager*  fEmManager;

protected:

  G4ElementData*               fElementData;
  G4VParticleChange*           pParticleChange;
  G4PhysicsTable*              xSectionTable;
  const std::vector<G4double>* theDensityFactor;
  const std::vector<G4int>*    theDensityIdx;
  size_t                       idxTable;
  G4bool                       lossFlucFlag;
  const static G4double        inveplus;       

  // ======== Cached values - may be state dependent ================

private:

  const G4MaterialCutsCouple* fCurrentCouple;
  const G4Element*            fCurrentElement;
  const G4Isotope*            fCurrentIsotope;
  G4VEmModel*                 fTripletModel;

  G4int                  nsec;
  std::vector<G4double>  xsec;

};

// ======== Run time inline methods ================

inline void G4VEmModel::SetCurrentCouple(const G4MaterialCutsCouple* p)
{
  fCurrentCouple = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4MaterialCutsCouple* G4VEmModel::CurrentCouple() const
{
  return fCurrentCouple;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetCurrentElement(const G4Element* elm)
{
  fCurrentElement = elm;
  fCurrentIsotope = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4Element* G4VEmModel::GetCurrentElement() const
{
  return fCurrentElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4Isotope* G4VEmModel::GetCurrentIsotope() const
{
  return fCurrentIsotope;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4VEmModel::MaxSecondaryKinEnergy(const G4DynamicParticle* dynPart)
{
  return MaxSecondaryEnergy(dynPart->GetParticleDefinition(),
                            dynPart->GetKineticEnergy());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmModel::ComputeDEDX(const G4MaterialCutsCouple* couple,
                                        const G4ParticleDefinition* part,
                                        G4double kinEnergy,
                                        G4double cutEnergy)
{
  SetCurrentCouple(couple);
  return ComputeDEDXPerVolume(couple->GetMaterial(),part,kinEnergy,cutEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmModel::CrossSection(const G4MaterialCutsCouple* couple,
                                         const G4ParticleDefinition* part,
                                         G4double kinEnergy,
                                         G4double cutEnergy,
                                         G4double maxEnergy)
{
  SetCurrentCouple(couple);
  return CrossSectionPerVolume(couple->GetMaterial(),part,kinEnergy,
                               cutEnergy,maxEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4VEmModel::ComputeMeanFreePath(const G4ParticleDefinition* part,
                                         G4double ekin,
                                         const G4Material* material,
                                         G4double emin,
                                         G4double emax)
{
  G4double cross = CrossSectionPerVolume(material,part,ekin,emin,emax);
  return cross > 0.0 ? 1./cross : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double
G4VEmModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition* part,
                                       const G4Element* elm,
                                       G4double kinEnergy,
                                       G4double cutEnergy,
                                       G4double maxEnergy)
{
  SetCurrentElement(elm);
  return ComputeCrossSectionPerAtom(part,kinEnergy,elm->GetZ(),elm->GetN(),
                                    cutEnergy,maxEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4Element*
G4VEmModel::SelectRandomAtom(const G4MaterialCutsCouple* couple,
                             const G4ParticleDefinition* part,
                             G4double kinEnergy,
                             G4double cutEnergy,
                             G4double maxEnergy)
{
  fCurrentCouple = couple;
  if(nSelectors > 0) {
    fCurrentElement = 
      ((*elmSelectors)[couple->GetIndex()])->SelectRandomAtom(kinEnergy);
  } else {
    fCurrentElement = SelectRandomAtom(couple->GetMaterial(),part,kinEnergy,
                                       cutEnergy,maxEnergy);
  }
  fCurrentIsotope = nullptr;
  return fCurrentElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEmModel::SelectRandomAtomNumber(const G4Material* mat)
{
  // this algorith assumes that cross section is proportional to
  // number electrons multiplied by number of atoms
  size_t nn = mat->GetNumberOfElements();
  const G4ElementVector* elmv = mat->GetElementVector();
  G4int Z = (*elmv)[0]->GetZasInt();
  if(1 < nn) {
    const G4double* at = mat->GetVecNbOfAtomsPerVolume();
    G4double tot = mat->GetTotNbOfAtomsPerVolume()*G4UniformRand();
    for( size_t i=0; i<nn; ++i) {
      tot -= at[i];
      if(tot <= 0.0) { 
	Z = (*elmv)[i]->GetZasInt();
	break; 
      }
    }
  }
  return Z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4VEmModel::SelectIsotopeNumber(const G4Element* elm)
{
  SetCurrentElement(elm);
  G4int N = G4lrint(elm->GetN());
  G4int ni = elm->GetNumberOfIsotopes();
  fCurrentIsotope = nullptr;
  if(ni > 0) {
    G4int idx = 0;
    if(ni > 1) {
      G4double* ab = elm->GetRelativeAbundanceVector();
      G4double x = G4UniformRand();
      for(; idx<ni; ++idx) {
        x -= ab[idx];
        if (x <= 0.0) { break; }
      }
      if(idx >= ni) { idx = ni - 1; }
    }
    fCurrentIsotope = elm->GetIsotope(idx);
    N = fCurrentIsotope->GetN();
  }
  return N;
}

// ======== Get/Set inline methods used at initialisation ================

inline G4VEmFluctuationModel* G4VEmModel::GetModelOfFluctuations()
{
  return flucModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmAngularDistribution* G4VEmModel::GetAngularDistribution()
{
  return anglModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetAngularDistribution(G4VEmAngularDistribution* p)
{
  if(p != anglModel) {
    delete anglModel;
    anglModel = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4VEmModel::GetTripletModel()
{
  return fTripletModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetTripletModel(G4VEmModel* p)
{
  if(p != fTripletModel) {
    delete fTripletModel;
    fTripletModel = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmModel::HighEnergyLimit() const
{
  return highLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmModel::LowEnergyLimit() const
{
  return lowLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmModel::HighEnergyActivationLimit() const
{
  return eMaxActive;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmModel::LowEnergyActivationLimit() const
{
  return eMinActive;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmModel::PolarAngleLimit() const
{
  return polarAngleLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmModel::SecondaryThreshold() const
{
  return secondaryThreshold;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEmModel::LPMFlag() const 
{
  return theLPMflag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEmModel::DeexcitationFlag() const 
{
  return flagDeexcitation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEmModel::ForceBuildTableFlag() const 
{
  return flagForceBuildTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEmModel::UseAngularGeneratorFlag() const
{
  return useAngularGenerator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetAngularGeneratorFlag(G4bool val)
{
  useAngularGenerator = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetFluctuationFlag(G4bool val)
{
  lossFlucFlag = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetMasterThread(G4bool val)
{
  isMaster = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEmModel::IsMaster() const
{
  return isMaster;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetHighEnergyLimit(G4double val)
{
  highLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetLowEnergyLimit(G4double val)
{
  lowLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetActivationHighEnergyLimit(G4double val)
{
  eMaxActive = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetActivationLowEnergyLimit(G4double val)
{
  eMinActive = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEmModel::IsActive(G4double kinEnergy)
{
  return (kinEnergy >= eMinActive && kinEnergy <= eMaxActive);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetPolarAngleLimit(G4double val)
{
  if(!isLocked) { polarAngleLimit = val; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetSecondaryThreshold(G4double val) 
{
  secondaryThreshold = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetLPMFlag(G4bool val) 
{
  theLPMflag = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetDeexcitationFlag(G4bool val) 
{
  flagDeexcitation = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetForceBuildTable(G4bool val)
{
  flagForceBuildTable = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4String& G4VEmModel::GetName() const 
{
  return name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4EmElementSelector*>* G4VEmModel::GetElementSelectors()
{
  return elmSelectors;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4VEmModel::SetElementSelectors(std::vector<G4EmElementSelector*>* p)
{
  elmSelectors = p;
  if(elmSelectors) { nSelectors = elmSelectors->size(); }
  localElmSelectors = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4ElementData* G4VEmModel::GetElementData()
{
  return fElementData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4PhysicsTable* G4VEmModel::GetCrossSectionTable()
{
  return xSectionTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4VEmModel::IsLocked() const
{
  return isLocked;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4VEmModel::SetLocked(G4bool val)
{
  isLocked = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif

