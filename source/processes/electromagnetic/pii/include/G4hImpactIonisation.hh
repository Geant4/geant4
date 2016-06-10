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
// ------------------------------------------------------------
// G4hImpactIonisation
//
// $Id: G4hImpactIonisation.hh 70904 2013-06-07 10:34:25Z gcosmo $
//
// Author: Maria Grazia Pia (MariaGrazia.Pia@ge.infn.it)
//
// 08 Sep 2008 - MGP - Created (initially based on G4hLowEnergyIonisation) 
//                     Added PIXE capabilities
//                     Partial clean-up of the implementation (more needed)
//                     Calculation of MicroscopicCrossSection delegated to specialised class 
//
// ------------------------------------------------------------
 
// Class Description:
// Impact Ionisation process of charged hadrons and ions
// Initially based on G4hLowEnergyIonisation, to be subject to redesign
// and further evolution of physics capabilities
//
// The physics model of G4hLowEnergyIonisation is described in 
// CERN-OPEN-99-121 and CERN-OPEN-99-300. 
//
// Documentation available in:
// M.G. Pia et al., PIXE Simulation With Geant4,
// IEEE Trans. Nucl. Sci., vol. 56, no. 6, pp. 3614-3649, Dec. 2009.

// ------------------------------------------------------------

#ifndef G4HIMPACTIONISATION
#define G4HIMPACTIONISATION 1

#include <map>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4hRDEnergyLoss.hh"
#include "G4DataVector.hh"
#include "G4AtomicDeexcitation.hh"
#include "G4PixeCrossSectionHandler.hh"

class G4VLowEnergyModel;
class G4VParticleChange;
class G4ParticleDefinition;
class G4PhysicsTable;
class G4MaterialCutsCouple;
class G4Track;
class G4Step;

class G4hImpactIonisation : public G4hRDEnergyLoss
{
public: // With description
  
  G4hImpactIonisation(const G4String& processName = "hImpactIoni"); 
  // The ionisation process for hadrons/ions to be include in the
  // UserPhysicsList

  ~G4hImpactIonisation();
  // Destructor
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  // True for all charged hadrons/ions
    
  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) ;
  // Build physics table during initialisation

  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   enum G4ForceCondition* condition );
  // Return MeanFreePath until delta-electron production
  
  void PrintInfoDefinition() const;
  // Print out of the class parameters

  void SetHighEnergyForProtonParametrisation(G4double energy) {protonHighEnergy = energy;} ;
  // Definition of the boundary proton energy. For higher energies
  // Bethe-Bloch formula is used, for lower energies a parametrisation
  // of the energy losses is performed. Default is 2 MeV.

  void SetLowEnergyForProtonParametrisation(G4double energy) {protonLowEnergy = energy;} ;
  // Set of the boundary proton energy. For lower energies
  // the Free Electron Gas model is used for the energy losses.
  // Default is 1 keV.

  void SetHighEnergyForAntiProtonParametrisation(G4double energy) {antiprotonHighEnergy = energy;} ;
  // Set of the boundary antiproton energy. For higher energies
  // Bethe-Bloch formula is used, for lower energies parametrisation
  // of the energy losses is performed. Default is 2 MeV.

  void SetLowEnergyForAntiProtonParametrisation(G4double energy) {antiprotonLowEnergy = energy;} ;
  // Set of the boundary antiproton energy. For lower energies
  // the Free Electron Gas model is used for the energy losses. 
  // Default is 1 keV.

  G4double GetContinuousStepLimit(const G4Track& track,
				  G4double previousStepSize,
				  G4double currentMinimumStep,
				  G4double& currentSafety); 
  // Calculation of the step limit due to ionisation losses

  void SetElectronicStoppingPowerModel(const G4ParticleDefinition* aParticle,
                                       const G4String& dedxTable);
  // This method defines the electron ionisation parametrisation method 
  // via the name of the table. Default is "ICRU_49p".

  void SetNuclearStoppingPowerModel(const G4String& dedxTable)
  {theNuclearTable = dedxTable; SetNuclearStoppingOn();};
  // This method defines the nuclear ionisation parametrisation method 
  // via the name of the table. Default is "ICRU_49".

  // ---- MGP ---- The following design of On/Off is nonsense; to be modified
  // in a following design iteration

  void SetNuclearStoppingOn() {nStopping = true;};
  // This method switch on calculation of the nuclear stopping power.
  
  void SetNuclearStoppingOff() {nStopping = false;};
  // This method switch off calculation of the nuclear stopping power.

  void SetBarkasOn() {theBarkas = true;};
  // This method switch on calculation of the Barkas and Bloch effects.

  void SetBarkasOff() {theBarkas = false;};
  // This method switch off calculation of the Barkas and Bloch effects.

  void SetPixe(const G4bool /* val */ ) {pixeIsActive = true;};
  // This method switches atomic relaxation on/off; currently always on

  G4VParticleChange* AlongStepDoIt(const G4Track& trackData ,
                                   const G4Step& stepData ) ;
  // Function to determine total energy deposition on the step

  G4VParticleChange* PostStepDoIt(const G4Track& track,
				  const G4Step& Step  ) ;
  // Simulation of delta-ray production.

  G4double ComputeDEDX(const G4ParticleDefinition* aParticle,
                       const G4MaterialCutsCouple* couple,
		       G4double kineticEnergy);
  // This method returns electronic dE/dx for protons or antiproton

  void SetCutForSecondaryPhotons(G4double cut);
  // Set threshold energy for fluorescence

  void SetCutForAugerElectrons(G4double cut);
  // Set threshold energy for Auger electron production

  void ActivateAugerElectronProduction(G4bool val);
  // Set Auger electron production flag on/off

  // Accessors to configure PIXE
  void SetPixeCrossSectionK(const G4String& name) { modelK = name; }
  void SetPixeCrossSectionL(const G4String& name) { modelL = name; }
  void SetPixeCrossSectionM(const G4String& name) { modelM = name; }
  void SetPixeProjectileMinEnergy(G4double energy) { eMinPixe = energy; }
  void SetPixeProjectileMaxEnergy(G4double energy) { eMaxPixe = energy; }

protected:

private:

  void InitializeMe();
  void InitializeParametrisation();
  void BuildLossTable(const G4ParticleDefinition& aParticleType);
  // void BuildDataForFluorescence(const G4ParticleDefinition& aParticleType);
  void BuildLambdaTable(const G4ParticleDefinition& aParticleType);
  void SetProtonElectronicStoppingPowerModel(const G4String& dedxTable)
  {protonTable = dedxTable ;};
  // This method defines the ionisation parametrisation method via its name

  void SetAntiProtonElectronicStoppingPowerModel(const G4String& dedxTable)
  {antiprotonTable = dedxTable;};

  G4double MicroscopicCrossSection(const G4ParticleDefinition& aParticleType,
				   G4double kineticEnergy,
				   G4double atomicNumber,
				   G4double deltaCutInEnergy) const;

  G4double GetConstraints(const G4DynamicParticle* particle,
                          const G4MaterialCutsCouple* couple);
  // Function to determine StepLimit

  G4double ProtonParametrisedDEDX(const G4MaterialCutsCouple* couple,
				  G4double kineticEnergy) const;

  G4double AntiProtonParametrisedDEDX(const G4MaterialCutsCouple* couple,
				      G4double kineticEnergy) const;

  G4double DeltaRaysEnergy(const G4MaterialCutsCouple* couple,
			   G4double kineticEnergy,
			   G4double particleMass) const;
  // This method returns average energy loss due to delta-rays emission with
  // energy higher than the cut energy for given material.

  G4double BarkasTerm(const G4Material* material,
		      G4double kineticEnergy) const;
  // Function to compute the Barkas term for protons

  G4double BlochTerm(const G4Material* material,
		     G4double kineticEnergy,
		     G4double cSquare) const;
  // Function to compute the Bloch term	for protons

  G4double ElectronicLossFluctuation(const G4DynamicParticle* particle,
                                     const G4MaterialCutsCouple* material,
				     G4double meanLoss,
				     G4double step) const;
  // Function to sample electronic losses

  // hide assignment operator
  G4hImpactIonisation & operator=(const G4hImpactIonisation &right);
  G4hImpactIonisation(const G4hImpactIonisation&);

private:
  //  private data members ...............................
  G4VLowEnergyModel* betheBlochModel;
  G4VLowEnergyModel* protonModel;
  G4VLowEnergyModel* antiprotonModel;
  G4VLowEnergyModel* theIonEffChargeModel;
  G4VLowEnergyModel* theNuclearStoppingModel;
  G4VLowEnergyModel* theIonChuFluctuationModel;
  G4VLowEnergyModel* theIonYangFluctuationModel;

  // std::map<G4int,G4double,std::less<G4int> > totalCrossSectionMap;

  // name of parametrisation table of electron stopping power
  G4String protonTable;
  G4String antiprotonTable;
  G4String theNuclearTable;

  // interval of parametrisation of electron stopping power
  G4double protonLowEnergy;
  G4double protonHighEnergy;
  G4double antiprotonLowEnergy;
  G4double antiprotonHighEnergy;

  // flag of parametrisation of nucleus stopping power
  G4bool nStopping;
  G4bool theBarkas;

  G4DataVector cutForDelta;
  G4DataVector cutForGamma;
  G4double minGammaEnergy;
  G4double minElectronEnergy;
  G4PhysicsTable* theMeanFreePathTable;

  const G4double paramStepLimit; // parameter limits the step at low energy

  G4double fdEdx;        // computed in GetContraints
  G4double fRangeNow ;   //
  G4double charge;       //
  G4double chargeSquare; //
  G4double initialMass;  // mass to calculate Lambda tables
  G4double fBarkas;

  G4PixeCrossSectionHandler* pixeCrossSectionHandler;
  G4AtomicDeexcitation atomicDeexcitation;
  G4String modelK;
  G4String modelL;
  G4String modelM;
  G4double eMinPixe;
  G4double eMaxPixe;
 
  G4bool pixeIsActive;

};


inline G4double G4hImpactIonisation::GetContinuousStepLimit(const G4Track& track,
							    G4double,
							    G4double currentMinimumStep,
							    G4double&)
{
  G4double step = GetConstraints(track.GetDynamicParticle(),track.GetMaterialCutsCouple()) ;

  // ---- MGP ---- The following line, taken as is from G4hLowEnergyIonisation,
  // is meaningless: currentMinimumStep is passed by value,
  // therefore any local modification to it has no effect

  if ((step > 0.) && (step < currentMinimumStep)) currentMinimumStep = step ;

  return step ;
}


inline G4bool G4hImpactIonisation::IsApplicable(const G4ParticleDefinition& particle)
{
  // ---- MGP ---- Better criterion for applicability to be defined;
  // now hard-coded particle mass > 0.1 * proton_mass

  return (particle.GetPDGCharge() != 0.0 && particle.GetPDGMass() > CLHEP::proton_mass_c2*0.1);
}

#endif








