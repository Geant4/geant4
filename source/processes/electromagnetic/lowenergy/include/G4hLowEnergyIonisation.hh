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

// ************************************************************
// 28 July 1999 V.Ivanchenko cleen up
// 17 August 1999 G.Mancinelli implemented ICRU parametrization (protons)  
// 20 August 1999 G.Mancinelli implemented ICRU parametrization (alpha)  
// 31 August 1999 V.Ivanchenko update and cleen up 
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
public: // Without description
  
  G4hLowEnergyIonisation(const G4String& processName = "hLowEIoni"); 
  
  ~G4hLowEnergyIonisation();
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
  
  void BuildLambdaTable(const G4ParticleDefinition& aParticleType);
  
  G4double GetMeanFreePath(
			   const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition ) ;
  
  G4VParticleChange *PostStepDoIt(const G4Track& track,
				  const G4Step& Step  ) ;                 

  void BuildLossTable(const G4ParticleDefinition& aParticleType);

  void PrintInfoDefinition();
  
protected:
  
  virtual G4double ComputeMicroscopicCrossSection(
						  const G4ParticleDefinition& aParticleType,
						  G4double KineticEnergy,
						  G4double AtomicNumber);
      
public: // With description
  
  void SetStoppingPowerTableName(const G4String& dedxTable);
  // This method defines the ionisation parametrisation method via its name 

  void SetNuclearStoppingOn();
  // This method switch on calculation of the nuclear stopping power.
  
  void SetNuclearStoppingOff();
  // This method switch off calculation of the nuclear stopping power.
  
  void SetAntiProtonStoppingOn();
  // This method switch on calculation of the Barkas Effect for antiproton
  
  void SetAntiProtonStoppingOff();
  // This method switch on calculation of the Barkas Effect for antiproton
  
  virtual G4double GetParametrisedLoss(G4Material* aMaterial,
  		        	       const G4double KinEnergy,
			               const G4double DeltaRayCutNow);
  // This method returns parametrised energy loss.

  G4double GetPreciseDEDX(G4Material* aMaterial,
  			  const G4double KinEnergy,
		          const G4ParticleDefinition* aParticleType);
  // This method returns electron ionisation energy loss for any energy.

  G4double GetNuclearDEDX(G4Material* aMaterial,
  			  const G4double KinEnergy,
		          const G4ParticleDefinition* aParticleType);
  // This method returns nuclear energy loss.
  
  G4double GetBetheBlochLoss(const G4Material* material, 
                             const G4double KinEnergy,
			     const G4double DeltaRayCutNow);
  // This method returns energy loss calculated via Bethe-Bloch formula.
  
  G4double GetFreeElectronGasLoss(G4double paramA, G4double KinEnergy);
  // This method returns energy loss parametrised in the free electron gas model.
  
  G4double GetUrbanModel(const G4Element* element, G4double KinEnergy);
  // This method returns energy loss parametrised as in the hIonisation class.
  
  G4double GetDeltaRaysEnergy(const G4Material* material, const G4double KinEnergy,
			      const G4double DeltaRayCutNow);
  // This method returns average energy loss due to delta-rays emission with 
  // energy higher than the cut energy for given material.
  
  G4int MolecIsInICRU_R49p(const G4Material*  material);
  // This method returns index of the material in the table of protons energy
  // loss in ICRU Report N49. If material is not in the table the method returns -1.

  G4int MolecIsInICRU_R49PowersHe(const G4Material*  material);
  // This method returns index of the material in the table of He energy loss
  // in ICRU Report N49. If material is not in the table the method returns -1.

  G4double MolecIsInZiegler1988(const G4Material*  material);
  // This method returns index of the material in the table of energy loss from
  // NIM B35 (1988) 215-228. If material is not in the table the method returns -1.

  G4double GetMolecICRU_R49Loss(const G4Material* material, const G4double KinEnergy, 
			        const G4double DeltaRayCutNow, const G4int molecIndex);
  // This method returns energy loss of protons in material from the table of ICRU 
  // Report N49.

  G4double GetChemicalFactor(const G4double ExpStopPower125, const G4double KinEnergy,
			     const G4double BraggStopPower125);
  // This method returns the value of "chemical factor" which allows to correct
  // energy losses calculated according to the Bragg's rule (NIM B35 (1988) 215-228).
  
  G4double GetStoppingPower1977H(G4int iz, G4double E);
  // This method returns protons electronic stopping power parametrised according to
  // H.H.Andersen & J.F.Ziegler, Hydrogen Stopping Powers and
  // Ranges in All Elements, Vol.3, Pergamon Press, 1977
  
  G4double GetStoppingPowerICRU_R49p(G4int iz, G4double E, G4String type);
  // This method returns protons electronic stopping power parametrised according to
  // ICRU Report N49, 1993.
  
  G4double GetStoppingPower1977He(G4int iz, G4double E);
  // This method returns He electronic stopping power parametrised according to
  // J.F.Ziegler, Helium Stopping Powers and
  // Ranges in All Elemental Matter, Vol.4, Pergamon Press, 1977
  
  G4double GetStoppingPowerICRU_R49He(G4int iz, G4double E);
  // This method returns He electronic stopping power parametrised according to
  // ICRU Report N49, 1993. J.F. Ziegler model.
  
  G4double GetStoppingPowerICRU_R49PowersHe(G4int iz, G4double E);
  // This method returns He electronic stopping power parametrised according to
  // J.F.Ziegler, Helium Stopping Powers and
  // Ranges in All Elemental Matter, Vol.4, Pergamon Press, 1977
  
  G4double GetStoppingPower1977n(G4double Z1, G4double Z2, 
				 G4double M1, G4double M2, G4double E);
  // This method returns nuclear stopping power parametrised according to
  // J.F.Ziegler, Helium Stopping Powers and
  // Ranges in All Elemental Matter, Vol.4, Pergamon Press, 1977
  
  G4double GetStoppingPower1985n(G4double Z1, G4double Z2, 
				 G4double M1, G4double M2, G4double E);
  // This method returns nuclear stopping power parametrised according to
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985
  
  G4double GetStoppingPowerMoliere(G4double Z1, G4double Z2, 
                                   G4double M1, G4double M2, G4double E);
  // This method returns nuclear stopping power parametrised according to
  // ICRU Report N49, 1993. Moliere model.
  
  G4double GetHeEffChargeSquare(const G4int iz, const G4double HeKinEnergy);
  // This method returns He effective charge square parametrised according to
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985

  G4double GetIonEffChargeSquare(const G4Material* material, const G4double KinEnergy,
                                 const G4double IonCharge);
  // This method returns ion effective charge square parametrised according to
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985

  G4double ComputeBarkasTerm(const G4Material* material, const G4double KinEnergy);
  // Function to compute the Barkas term						  

  G4double GetConstraints(const G4DynamicParticle *aParticle,
                          G4Material *aMaterial);
  // Function to determine StepLimit
                                       
  G4VParticleChange* AlongStepDoIt(const G4Track& trackData , 
                                   const G4Step& stepData );
  // Function to determine total energy deposition on the step

private:
  
  // hide assignment operator 
  G4hLowEnergyIonisation & operator=(const G4hLowEnergyIonisation &right);
  G4hLowEnergyIonisation(const G4hLowEnergyIonisation&);
  
private:
  //  private data members ...............................

protected:
  //  protected data members ...............................
  
  G4PhysicsTable* theMeanFreePathTable;
  
  // interval of parametrisation of electron stopping power 
  G4double ParamLowEnergy;
  G4double ParamHighEnergy;
  // name of parametrisation table of electron stopping power
  G4String DEDXtable;
  // flag of parametrisation of nucleus stopping power
  G4bool nStopping;
  G4bool pbarStop;
  
  // constants needed for the energy loss calculation
  
  const G4double twoln10;
  const G4double Factor;
  const G4double bg2lim;
  const G4double taulim;          // energy to start to switch off shell corrections
  G4double RateMass;              // m_e/M
  G4double MassRatio;             // m_p/M
  G4ParticleDefinition* theParticle;
  
  // particles , cuts in kinetic energy ........
  const G4Electron* theElectron;
  const G4Proton* theProton;
  const G4AntiProton* theAntiProton;
  
  const G4double* DeltaCutInKineticEnergy ; 
  
  G4double DeltaCutInKineticEnergyNow ;
  
  G4double ProtonMassAMU;
  G4double HeMassAMU;
  G4double ZieglerFactor; // Factor to convert the Stopping Power 
  // unit [ev/(10^15 atoms/cm^2]
  // into the Geant4 dE/dx unit
    
   static G4double LowerBoundLambda ; // bining for lambda table
   static G4double UpperBoundLambda ;
   static G4int    NbinLambda ;

   G4double LowestKineticEnergy,HighestKineticEnergy ;
   G4int    TotBin ;

  public:

    static void SetLowerBoundLambda(G4double val) {LowerBoundLambda = val;};
    static void SetUpperBoundLambda(G4double val) {UpperBoundLambda = val;};
    static void SetNbinLambda(G4int n) {NbinLambda = n;};
    static G4double GetLowerBoundLambda() { return LowerBoundLambda;};
    static G4double GetUpperBoundLambda() { return UpperBoundLambda;};
    static G4int GetNbinLambda() {return NbinLambda;};

};

#include "G4hLowEnergyIonisation.icc"

#endif
 







