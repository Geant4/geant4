//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hLowEnergyIonisationMA physics process -----
//                by Vladimir Ivanchenko, 14 July 1999 
//                was made on the base of G4hIonisation class
//                developed by Laszlo Urban
// ************************************************************

// ************************************************************
// 28 July   1999 V.Ivanchenko cleen up
// 17 August 1999 G.Mancinelli implemented ICRU parametrization (protons)  
// 20 August 1999 G.Mancinelli implemented ICRU parametrization (alpha)  
// 31 August 1999 V.Ivanchenko update and cleen up 
// 23 May    2000    MG Pia  Clean up for QAO model
// 25 July   2000 V.Ivanchenko New design iteration
// 09 August 2000 V.Ivanchenko Add GetContinuousStepLimit
// 17 August 2000 V.Ivanchenko Add IonFluctuationModel
// 23 Oct    2000 V.Ivanchenko Renew comments
// 30 Oct    2001 V.Ivanchenko Add minGammaEnergy and minElectronEnergy
// 07 Dec    2001 V.Ivanchenko Add SetFluorescence method
// 26 Feb    2002 V.Ivanchenko Add initialMass for GenericIons
// 21 Jan    2003 V.Ivanchenko Cut per region
// 17 April  2004 V.Ivanchenko migrate to model design
// ------------------------------------------------------------

// Class Description:
// Ionisation process of charged hadrons and ions, including low energy
// extensions
// The physics model is described in CERN-OPEN-99-121 and CERN-OPEN-99-300.
// The user may select parametrisation tables for electronic
// stopping powers and nuclear stopping powers
// The list of available tables:
// Electronic stopping powers: "ICRU_49p" (default), "ICRU_49He",
//                             "Ziegler1977p", "Ziegler1985p",
//                             "Ziegler1977He"
// Nuclear stopping powers:    "ICRU_49" (default), "Ziegler1977",
//                             "Ziegler1985"
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// and in the Physics Reference Manual

// ------------------------------------------------------------

#ifndef G4hLowEnergyIonisationMA_h
#define G4hLowEnergyIonisationMA_h 1

#include "G4VEnergyLossProcess.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4DataVector.hh"
#include "G4AtomicDeexcitation.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ionEffectiveCharge.hh"

class G4ShellVacancy;
class G4VhShellCrossSection;
class G4VEMDataSet;
class G4Region;
class G4Step;

class G4hLowEnergyIonisationMA : public G4VEnergyLossProcess
{
public: // With description

  G4hLowEnergyIonisationMA(const G4String& processName = "hLowEIoni");
  // The ionisation process for hadrons/ions to be include in the
  // UserPhysicsList

  ~G4hLowEnergyIonisationMA();
  // Destructor

  G4bool IsApplicable(const G4ParticleDefinition&);
  // True for all charged hadrons/ions

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
                                    const G4Material*, G4double cut);

  virtual std::vector<G4Track*>* SecondariesAlongStep(
                             const G4Step&,
                                   G4double&,
                                   G4double&,
                                   G4double&);

  virtual void SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double&,
                                   G4double&);

  virtual void PrintInfoDefinition();
  // Print out of the class parameters

  void SetHighEnergyForParametrisation(G4double energy)
                             {highEnergy = energy;} ;
  // Set of the boundary antiproton energy. For higher energies
  // Bethe-Bloch formula is used, for lower energies parametrisation
  // of the energy losses is performed. Default is 10 MeV.

  void SetElectronicStoppingPowerModel(const G4String& dedxTable)
                             {theTable = dedxTable;};
  // This method defines the electron ionisation parametrisation method
  // via the name of the table. Default is "ICRU_49p".

  void SetBarkasOn() {theBarkas = true;};
  // This method switch on calculation of the Barkas and Bloch effects.

  void SetBarkasOff() {theBarkas = false;};
  // This method switch off calculation of the Barkas and Bloch effects.

  virtual void ActivateFluorescence(G4bool val, const G4Region* r=0);
  void SetFluorescence(G4bool val) {ActivateFluorescence(val);};
  // This method switch on/off simulation of the fluorescence of the media.

  virtual void ActivateAugerElectronProduction(G4bool val, const G4Region* r=0);
  // Set Auger electron production flag on/off

  void SetCutForSecondaryPhotons(G4double cut) {minGammaEnergy = cut;};
  // Set threshold energy for fluorescence

  void SetCutForAugerElectrons(G4double cut) {minElectronEnergy = cut;};
  // Set threshold energy for Auger electron production

protected:

  const G4ParticleDefinition* DefineBaseParticle(const G4ParticleDefinition* p);

  virtual G4double GetMeanFreePath(const G4Track& track,
                                         G4double previousStepSize,
                                         G4ForceCondition* condition);

  virtual G4double GetContinuousStepLimit(const G4Track& track,
                                                G4double previousStepSize,
                                                G4double currentMinimumStep,
                                                G4double& currentSafety);

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dynParticle);

private:

  void InitialiseProcess();

  void BuildDataForFluorescence();

  G4double BarkasTerm(const G4Material* material,
                            G4double kineticEnergy) const;
  // Function to compute the Barkas term for protons

  G4double BlochTerm(const G4Material* material,
                           G4double kineticEnergy,
                           G4double cSquare) const;
  // Function to compute the Bloch term	for protons


  std::vector<G4DynamicParticle*>* DeexciteAtom(
                         const G4MaterialCutsCouple* couple,
 			       G4double incidentEnergy,
			       G4double tmax,
			       G4double hMass,
			       G4double eLoss);

  G4int SelectRandomAtom(const G4MaterialCutsCouple* couple,
                               G4double kineticEnergy) const;

  // hide assignment operator
  G4hLowEnergyIonisationMA & operator=(const G4hLowEnergyIonisationMA &right);
  G4hLowEnergyIonisationMA(const G4hLowEnergyIonisationMA&);

private:

  // name of parametrisation table of electron stopping power
  G4String                     theTable;

  // interval of parametrisation of electron stopping power
  G4double                     highEnergy;

  G4DataVector                 cutForDelta;
  G4DataVector                 cutForGamma;

  G4ionEffectiveCharge         effCharge;

  G4VEmFluctuationModel*       flucModel;
  G4AtomicDeexcitation         deexcitationManager;
  G4ShellVacancy*              shellVacancy;
  G4VhShellCrossSection*       shellCS;
  std::vector<G4VEMDataSet*>   zFluoDataVector;
  std::vector<const G4Region*> regionsWithFluo;

  // cash
  const G4Material*            theMaterial;
  const G4ParticleDefinition*  currentParticle;
  const G4ParticleDefinition*  theParticle;
  const G4ParticleDefinition*  theBaseParticle;
  G4double                     minGammaEnergy;
  G4double                     minElectronEnergy;

  G4double                     chargeCorrection;
  G4double                     chargeLowLimit;
  G4double                     energyLowLimit;
  G4double                     paramStepLimit;

  size_t                       fluobins;
  G4bool                       theBarkas;
  G4bool                       theFluo;
  G4bool                       fluoIsInitialised;
  G4bool                       isInitialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4hLowEnergyIonisationMA::IsApplicable(
                                const G4ParticleDefinition& particle)
{
  return(particle.GetPDGCharge() != 0.0
      && particle.GetPDGMass() > proton_mass_c2*0.1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4hLowEnergyIonisationMA::MinPrimaryEnergy(
          const G4ParticleDefinition*, const G4Material*, G4double cut)
{
  G4double x = 0.5*cut/electron_mass_c2;
  G4double g = sqrt(1. + x);
  return proton_mass_c2*(g - 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4hLowEnergyIonisationMA::MaxSecondaryEnergy(
          const G4DynamicParticle* dynParticle)
{
  G4double mass  = dynParticle->GetMass();
  G4double gamma = dynParticle->GetKineticEnergy()/mass + 1.0;
  G4double r     = electron_mass_c2/mass;
  G4double tmax  = 2.*electron_mass_c2*(gamma*gamma - 1.)/(1. + 2.*gamma*r + r*r);
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4hLowEnergyIonisationMA::GetMeanFreePath(
         const G4Track& track, G4double, G4ForceCondition* cond)
{
  G4double mRatio = proton_mass_c2/track.GetDynamicParticle()->GetMass();
  currentParticle = track.GetDefinition();
  theMaterial     = track.GetMaterial();
  G4double q_2    = effCharge.EffectiveChargeSquareRatio(currentParticle,theMaterial,
                                                         track.GetKineticEnergy());

  SetMassRatio(mRatio);
  SetReduceFactor(1.0/(q_2*mRatio));
  SetChargeSquare(q_2);
  SetChargeSquareRatio(q_2);
  return G4VEnergyLossProcess::GetMeanFreePath(track, 0.0, cond);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

