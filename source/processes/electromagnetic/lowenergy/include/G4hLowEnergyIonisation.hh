// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hLowEnergyIonisation physics process -----
//                by Vladimir Ivanchenko, 14 July 1999 
//                was made on the base of G4hIonisation class
//                developed by Laszlo Urban
// ************************************************************
// Class Description:
// G4hLowEnergyIonisation class is the extention of the ionisation 
// process for the slow charged hadrons. The physics model is
// described in CERN-OPEN-99-121. User have a possibility to define
// a parametrisation table via its name. 
// Class Description - End
//
// ************************************************************
// 23 May 2000    MG Pia  Clean up for QAO model 
// 28 July 1999 V.Ivanchenko cleen up
// 17 August 1999 G.Mancinelli implemented ICRU parametrization (protons)  
// 20 August 1999 G.Mancinelli implemented ICRU parametrization (alpha)  
// 31 August 1999 V.Ivanchenko update and cleen up 
// 25 July   2000 V.Ivanchenko New design iteration
// 09 July   2000 V.Ivanchenko Add GetContinuousStepLimit
// ------------------------------------------------------------
 
#ifndef G4hLowEnergyIonisation_h
#define G4hLowEnergyIonisation_h 1
 
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

class G4hLowEnergyIonisation : public G4hLowEnergyLoss
{
public: // Without description
  
  G4hLowEnergyIonisation(const G4String& processName = "hLowEIoni"); 

  ~G4hLowEnergyIonisation();
  
  G4bool IsApplicable(const G4ParticleDefinition&) const;
    
  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) ;

  void BuildLossTable(const G4ParticleDefinition& aParticleType) ;

  void BuildLambdaTable(const G4ParticleDefinition& aParticleType) ;
  
  G4double GetMeanFreePath(const G4Track& track,
			         G4double previousStepSize,
			    enum G4ForceCondition* condition );
  
  void PrintInfoDefinition() const;

  void SetHighEnergyForProtonParametrisation(G4double energy) 
                             {protonHighEnergy = energy;} ;

  void SetLowEnergyForProtonParametrisation(G4double energy) 
                             {protonLowEnergy = energy;} ;

  void SetHighEnergyForAntiProtonParametrisation(G4double energy) 
                             {antiProtonHighEnergy = energy;} ;

  void SetLowEnergyForAntiProtonParametrisation(G4double energy) 
                              {antiProtonLowEnergy = energy;} ;

  G4double GetContinuousStepLimit(const G4Track& track,
                                        G4double previousStepSize,
                                        G4double currentMinimumStep,
                                        G4double& currentSafety); 

public: // With description
  void SetElectronicStoppingPowerModel(const G4ParticleDefinition* aParticle,
                                       const G4String& dedxTable);
  // This method defines the ionisation parametrisation method via its name 

  void SetNuclearStoppingPowerModel(const G4String& dedxTable)
                 {theNuclearTable = dedxTable; SetNuclearStoppingOn();};
  // This method defines the ionisation parametrisation method via its name 

  void SetNuclearStoppingOn() {nStopping = true;};
  // This method switch on calculation of the nuclear stopping power.
  
  void SetNuclearStoppingOff() {nStopping = false;};
  // This method switch off calculation of the nuclear stopping power.
  
  void SetBarkasOn() {theBarkas = true;};
  // This method switch on calculation of the Barkas Effect for antiproton
  
  void SetBarkasOff() {theBarkas = false;};
  // This method switch on calculation of the Barkas Effect for antiproton
                                       
  G4VParticleChange* AlongStepDoIt(const G4Track& trackData , 
                                   const G4Step& stepData ) ;
  // Function to determine total energy deposition on the step

  G4VParticleChange* PostStepDoIt(const G4Track& track,
				  const G4Step& Step  ) ;                 
  // Simulation of delta rays production
    
  G4double ComputeDEDX(const G4ParticleDefinition* aParticle,
                       const G4Material* material,
                             G4double kineticEnergy);
  // This method returns total energy loss for a static particle.

protected:

private:

  void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);

  void InicialiseParametrisation();

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
  // Function to compute the Barkas term						  
  G4double BlochTerm(const G4Material* material,
                           G4double kineticEnergy,
                           G4double particleMass, 
                           G4double chargeSquare) const; 
  // Function to compute the Bloch term	

  G4double ElectronicLossFluctuation(const G4DynamicParticle* particle,
                                     const G4Material* material,
                                           G4double chargeSquare,
                                           G4double meanLoss,
                                           G4double step) const;
  // Function to sample electronic losses
		    
  // hide assignment operator 
  G4hLowEnergyIonisation & operator=(const G4hLowEnergyIonisation &right);
  G4hLowEnergyIonisation(const G4hLowEnergyIonisation&);
  
private:
  //  private data members ...............................
  G4VLowEnergyModel* theBetheBlochModel;
  G4VLowEnergyModel* theProtonModel;
  G4VLowEnergyModel* theAntiProtonModel;
  G4VLowEnergyModel* theIonEffChargeModel;
  G4VLowEnergyModel* theNuclearStoppingModel;

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

  G4double* deltaCutInKineticEnergy;
  G4PhysicsTable* theMeanFreePathTable;
  
  // constants needed for the energy loss calculation  
  const G4double factor;
  const G4double protonMass;             
  
  // particles 
  const G4Electron* theElectron;
  const G4Proton* theProton;
  const G4AntiProton* theAntiProton;

  G4double fdEdx;      // computed in GetContraints
  G4double fRangeNow ; // computed in GetContraints
 
protected:

private:

};

#include "G4hLowEnergyIonisation.icc"

#endif
 







