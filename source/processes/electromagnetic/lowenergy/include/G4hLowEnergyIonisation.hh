// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
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
// It is the extention of the ionisation process for the slow 
// charged hadrons.
// ************************************************************
// 28 July 1999 V.Ivanchenko cleen up
// 17 August 1999 G.Mancinelli implemented ICRU parametrization (protons)  
// 20 August 1999 G.Mancinelli implemented ICRU parametrization (alpha)  
// ------------------------------------------------------------
 
#ifndef G4hLowEnergyIonisation_h
#define G4hLowEnergyIonisation_h 1
 
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4hEnergyLoss.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"

class G4hLowEnergyIonisation : public G4hEnergyLoss
{
public:
  
  G4hLowEnergyIonisation(const G4String& processName = "hLowEIoni"); 
  
  ~G4hLowEnergyIonisation();
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);
  
  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
  
  void BuildLambdaTable(const G4ParticleDefinition& aParticleType);
  
  G4double GetMeanFreePath(
			   const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition ) ;
  
  G4VParticleChange *PostStepDoIt(const G4Track& track,
				  const G4Step& Step  ) ;                 
protected:
  
  virtual G4double ComputeMicroscopicCrossSection(
						  const G4ParticleDefinition& aParticleType,
						  G4double KineticEnergy,
						  G4double AtomicNumber);
  
public:
  
  void BuildLossTableZiegler(const G4ParticleDefinition& aParticleType);
  
  void BuildLossTableICRU(const G4ParticleDefinition& aParticleType);
  
  void SetStoppingPowerTableName(const G4String& dedxTable);
  
  void SetNuclearStoppingOn();
  
  void SetNuclearStoppingOff();
  
  G4double GetZieglerLoss(const G4Material* material, const G4double KinEnergy, 
			  const G4double DeltaRayCutNow);
  
  G4double GetICRULoss(const G4Material* material, const G4double KinEnergy, 
		       const G4double DeltaRayCutNow);
  
  G4double GetBetheBlochLoss(const G4Material* material, const G4double KinEnergy,
			     const G4double DeltaRayCutNow);
  
  G4double GetFreeElectronGasLoss(G4double paramA, G4double KinEnergy);
  
  G4double GetStoppingPower1977H(G4int iz, G4double E);
  
  G4double GetStoppingPowerICRU_R49p(G4int iz, G4double E, G4String type);
  
  G4double GetStoppingPower1977He(G4int iz, G4double E);
  
  G4double GetStoppingPowerICRU_R49He(G4int iz, G4double E);
  
  G4double GetStoppingPowerICRU_R49PowersHe(G4int iz, G4double E);
  
  G4double GetStoppingPower1977n(G4double Z1, G4double Z2, 
				 G4double M1, G4double M2, G4double E);
  
  G4double GetStoppingPower1985n(G4double Z1, G4double Z2, 
				 G4double M1, G4double M2, G4double E);
  
  G4double GetStoppingPowerMoliere(G4double Z1, G4double Z2, 
                                   G4double M1, G4double M2, G4double E);
  
  G4double GetUrbanModel(const G4Element* element, G4double KinEnergy);
  
  G4double GetDeltaRaysEnergy(const G4Material* material, const G4double KinEnergy,
			      const G4double DeltaRayCutNow);
  
  G4double GetChemicalFactor(const G4Material* material, const G4double KinEnergy,
			     const G4double DeltaRayCutNow);
  
  G4int MolecIsInICRUTable(const G4Material*  material);
  void PrintInfoDefinition();
  
private:
  
  // hide assignment operator 
  G4hLowEnergyIonisation & operator=(const G4hLowEnergyIonisation &right);
  G4hLowEnergyIonisation(const G4hLowEnergyIonisation&);
  
private:
  //  private data members ...............................
  
  G4PhysicsTable* theMeanFreePathTable;
  
  // interval of parametrisation of electron stopping power 
  G4double ParamLowEnergy;
  G4double ParamHighEnergy;
  // name of parametrisation table of electron stopping power
  G4String DEDXtable;
  // name of parametrisation type (Ziegler, ICRU)
  G4String Parametrization;
  // flag of parametrisation of nucleus stopping power
  G4bool nStopping;
  
  // constants needed for the energy loss calculation
  
  const G4double twoln10;
  const G4double Factor;
  const G4double bg2lim;
  const G4double taulim;          // energy to start to switch off shell corrections
  const G4double RateMass;
  
  // particles , cuts in kinetic energy ........
  const G4Electron* theElectron;
  const G4Proton* theProton;
  const G4AntiProton* theAntiProton;
  
  const G4double* DeltaCutInKineticEnergy ; 
  
  G4double DeltaCutInKineticEnergyNow ;
  
  G4double ProtonMassAMU;
  G4double ZieglerFactor; // Factor to convert the Stopping Power 
  // unit [ev/(10^15 atoms/cm^2]
  // into the Geant4 dE/dx unit
  
  
};

#include "G4hLowEnergyIonisation.icc"

#endif
 







