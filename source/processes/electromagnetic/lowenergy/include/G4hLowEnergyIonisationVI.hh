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
//      ---------- G4hLowEnergyIonisation physics process -----
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

#ifndef G4hLowEnergyIonisationVI_h
#define G4hLowEnergyIonisationVI_h 1
 
#include "G4hLowEnergyLoss.hh"
#include "G4VLowEnergyModel.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4hNuclearStoppingModel.hh"
#include "G4hBetheBlochModel.hh"
#include "G4hParametrisedLossModel.hh"
#include "G4QAOLowEnergyLoss.hh"
#include "G4hIonEffChargeSquare.hh"
#include "G4IonChuFluctuationModel.hh"
#include "G4IonYangFluctuationModel.hh"
#include "G4AtomicDeexcitation.hh"

class G4VEMDataSet;
class G4ShellVacancy;
class G4VhShellCrossSection;

class G4hLowEnergyIonisationVI : public G4hLowEnergyLoss
{
public: // With description
  
  G4hLowEnergyIonisationVI(const G4String& processName = "hLowEIoni"); 
  // The ionisation process for hadrons/ions to be include in the
  // UserPhysicsList

  ~G4hLowEnergyIonisationVI();
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

  void SetHighEnergyForProtonParametrisation(G4double energy) 
                             {protonHighEnergy = energy;} ;
  // Definition of the boundary proton energy. For higher energies
  // Bethe-Bloch formula is used, for lower energies a parametrisation
  // of the energy losses is performed. Default is 2 MeV.

  void SetLowEnergyForProtonParametrisation(G4double energy) 
                             {protonLowEnergy = energy;} ;
  // Set of the boundary proton energy. For lower energies
  // the Free Electron Gas model is used for the energy losses.
  // Default is 1 keV.

  void SetHighEnergyForAntiProtonParametrisation(G4double energy) 
                             {antiProtonHighEnergy = energy;} ;
  // Set of the boundary antiproton energy. For higher energies
  // Bethe-Bloch formula is used, for lower energies parametrisation
  // of the energy losses is performed. Default is 2 MeV.

  void SetLowEnergyForAntiProtonParametrisation(G4double energy) 
                              {antiProtonLowEnergy = energy;} ;
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

  void SetNuclearStoppingOn() {nStopping = true;};
  // This method switch on calculation of the nuclear stopping power.
  
  void SetNuclearStoppingOff() {nStopping = false;};
  // This method switch off calculation of the nuclear stopping power.
  
  void SetBarkasOn() {theBarkas = true;};
  // This method switch on calculation of the Barkas and Bloch effects. 
  
  void SetBarkasOff() {theBarkas = false;};
  // This method switch off calculation of the Barkas and Bloch effects.
                                       
  G4VParticleChange* AlongStepDoIt(const G4Track& trackData , 
                                   const G4Step& stepData ) ;
  // Function to determine total energy deposition on the step

  G4VParticleChange* PostStepDoIt(const G4Track& track,
				  const G4Step& Step  ) ;                 
  // Simulation of delta rays production.
    
  G4double ComputeDEDX(const G4ParticleDefinition* aParticle,
                       const G4Material* material,
                             G4double kineticEnergy);
  // This method returns electronic dE/dx for protons or antiproton.

protected:

private:

  void InitializeMe();

  void InitializeParametrisation();

  void BuildLossTable(const G4ParticleDefinition& aParticleType) ;

  void BuildLambdaTable(const G4ParticleDefinition& aParticleType) ;
  
  void SetProtonElectronicStoppingPowerModel(const G4String& dedxTable) 
                              {theProtonTable = dedxTable ;};
  // This method defines the ionisation parametrisation method via its name 

  void SetAntiProtonElectronicStoppingPowerModel(const G4String& dedxTable) 
                              {theAntiProtonTable = dedxTable ;};
  
  G4double ComputeMicroscopicCrossSection(
                  const G4ParticleDefinition& aParticleType,
	  	        G4double kineticEnergy,
			G4double atomicNumber,
                        G4double deltaCutInEnergy) const;

  G4double GetConstraints(const G4DynamicParticle* particle,
                          const G4Material* material);
  // Function to determine StepLimit

  G4double ProtonParametrisedDEDX(const G4Material* material, 
                                        G4double kineticEnergy) const;

  G4double AntiProtonParametrisedDEDX(const G4Material* material, 
                                            G4double kineticEnergy) const;
    
  G4double DeltaRaysEnergy(const G4Material* material, 
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
                                     const G4Material* material,
                                           G4double meanLoss,
                                           G4double step) const;
  // Function to sample electronic losses

  G4std::vector<G4DynamicParticle*>* DeexciteAtom(const G4Material* material,
						  G4double incidentEnergy,
						  G4double eLoss,
						  G4double hMass,
						  G4double hMomentum);

  G4int SelectRandomAtom(const G4Material* material, 
                               G4double kineticEnergy) const;
		    
  // hide assignment operator 
  G4hLowEnergyIonisationVI & operator=(const G4hLowEnergyIonisationVI &right);
  G4hLowEnergyIonisationVI(const G4hLowEnergyIonisationVI&);
  
private:
  //  private data members ...............................
  G4VLowEnergyModel* theBetheBlochModel;
  G4VLowEnergyModel* theProtonModel;
  G4VLowEnergyModel* theAntiProtonModel;
  G4VLowEnergyModel* theIonEffChargeModel;
  G4VLowEnergyModel* theNuclearStoppingModel;
  G4VLowEnergyModel* theIonChuFluctuationModel;
  G4VLowEnergyModel* theIonYangFluctuationModel;

  // name of parametrisation table of electron stopping power
  G4String theProtonTable;
  G4String theAntiProtonTable;
  G4String theNuclearTable;

  // interval of parametrisation of electron stopping power 
  G4double protonLowEnergy;
  G4double protonHighEnergy;
  G4double antiProtonLowEnergy;
  G4double antiProtonHighEnergy;

  // flag of parametrisation of nucleus stopping power
  G4bool nStopping;
  G4bool theBarkas;

  G4DataVector cutForDelta;
  G4double* gammaCutInEnergy;
  G4PhysicsTable* theMeanFreePathTable;
  
  const G4double paramStepLimit; // parameter limits the step at low energy
  
  G4double fdEdx;        // computed in GetContraints
  G4double fRangeNow ;   //         
  G4double charge;       //
  G4double chargeSquare; //
 
  G4AtomicDeexcitation deexcitationManager;
  G4ShellVacancy* shellVacancy;
  G4VhShellCrossSection* shellCS;
  G4std::vector<G4VEMDataSet*> zFluoDataVector;

};

#include "G4hLowEnergyIonisationVI.icc"

#endif
 







