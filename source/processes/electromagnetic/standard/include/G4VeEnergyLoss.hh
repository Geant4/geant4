// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VeEnergyLoss.hh,v 1.5 2000-10-30 06:56:19 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ---------------------------------------------------------------
// 18/11/98 It is a modified version of G4VeEnergyLoss: continuous energy loss
//          with generation of subcutoff delta rays, L. Urban
// 02/02/99 new data members ,  L.Urban 
// 10/02/00 modifications , new e.m. structure, L.Urban
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
#ifndef G4VeEnergyLoss_h
#define G4VeEnergyLoss_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VEnergyLoss.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4VParticleChange.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4EnergyLossTables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
// Class description:
//
// This class is the implementation of the unified Energy Loss process
// with generation of subcutoff secondaries , see description of 
// *****************************************
//  the AlongStepDoIt method.
// It calculates the continuous energy loss for e+/e-.
// The following processes give contributions to the continuous
// energy loss (by default) :
//  ---	 ionisation (= cont.ion.loss + delta ray production)
//  --- bremsstrahlung (= cont.loss due to soft brems+discrete bremsstrahlung)
//   more can be added	 ..........
// This class creates static dE/dx and range tables for e+ and e-,
// which tables can be used by other processes , too.
// G4VeEnergyLoss is the base class for the processes giving contribution
// to the (continuous) energy loss of e+/e- .
//
// Class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
class G4VeEnergyLoss : public G4VEnergyLoss
 
{
  public:
 
    G4VeEnergyLoss(const G4String& );

    virtual ~G4VeEnergyLoss();

  public: // With description

    G4bool IsApplicable(const G4ParticleDefinition&);
    //	true for e+/e- , false otherwise

    void BuildDEDXTable(const G4ParticleDefinition& aParticleType);
    //	It builds dE/dx and range tables for aParticleType and
    //	for every material contained in the materialtable.

    G4double GetContinuousStepLimit(const G4Track& track,
                                    G4double previousStepSize,
                                    G4double currentMinimumStep,
                                    G4double& currentSafety);
    // Computes the steplimit due to the energy loss process.

    G4VParticleChange* AlongStepDoIt(const G4Track& track,
                                     const G4Step& Step) ;
    // --- Performs the computation of the mean energy loss
    //     after the step .
    // --- Compares the value of the safety to the value
    //     of the cut in range.
    // --- If safety < cut in range , generates secondaries
    //     with kinetic energy in the interval [Tsafety,Tcut],
    //     where  Tsafety corresponds to range safety,
    //            Tcut is the cut in energy.
    // --- Updates the value of the mean energy loss and
    //     the direction of the primary
    // --- Computes the energy loss fluctuation.   

    virtual G4double GetMeanFreePath(const G4Track& track,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition) = 0;
    // Virtual function to be overridden in the derived classes
    // ( ionisation and bremsstrahlung) .

    virtual G4VParticleChange* PostStepDoIt(const G4Track& track,
                                            const G4Step& Step) = 0;
    // Virtual function to be overridden in the derived classes
    // ( ionisation and bremsstrahlung) .
                                            
                                            
  private:


    G4double GetConstraints(const G4DynamicParticle* aParticle,
                            G4Material* aMaterial);
                                                                  
    // hide  assignment operator
    G4VeEnergyLoss (G4VeEnergyLoss &); 
    G4VeEnergyLoss & operator=(const G4VeEnergyLoss &right);

  protected:

    G4PhysicsTable* theLossTable;
     
    G4double MinKineticEnergy ;           //

    G4double Charge,lastCharge ;


  private:

    G4PhysicsTable* theDEDXTable;

    G4int            CounterOfProcess;
    G4PhysicsTable** RecorderOfProcess;
                                            
    G4double fdEdx;                       // computed in GetConstraints
    G4double fRangeNow;                   // computed in GetConstraints

    G4double linLossLimit ;               //

    G4double c1N,c2N ;                   // coeffs to compute nb of deltas
    G4int Ndeltamax ;                    // upper limit for nb of subcutoff
                                         // delta rays in one step
    
 //  
 // static part of the class
 //
  public:  // With description

    static void  SetNbOfProcesses(G4int nb) {NbOfProcesses=nb;};
    // Sets number of processes giving contribution to the energy loss

    static void  PlusNbOfProcesses()        {NbOfProcesses++ ;};
    // Increases number of processes giving contribution to the energy loss

    static void  MinusNbOfProcesses()       {NbOfProcesses-- ;};
    // Decreases number of processes giving contribution to the energy loss

    static G4int GetNbOfProcesses()         {return NbOfProcesses;};
    // Gets number of processes giving contribution to the energy loss
    // ( default value = 2)

    static void SetLowerBoundEloss(G4double val) {LowerBoundEloss=val;};
    static void SetUpperBoundEloss(G4double val) {UpperBoundEloss=val;};
    static void SetNbinEloss(G4int nb)           {NbinEloss=nb;};

    static G4double GetLowerBoundEloss() {return LowerBoundEloss;};
    static G4double GetUpperBoundEloss() {return UpperBoundEloss;};
    static G4int    GetNbinEloss()       {return NbinEloss;};

 
 protected:
 
    //basic DEDX and Range tables
    static G4PhysicsTable* theDEDXElectronTable ;
    static G4PhysicsTable* theDEDXPositronTable ;
    static G4PhysicsTable* theRangeElectronTable ;
    static G4PhysicsTable* theRangePositronTable ;

    //inverse tables of the range tables
    static G4PhysicsTable* theInverseRangeElectronTable;
    static G4PhysicsTable* theInverseRangePositronTable;
   
    // lab and proper time tables
    static G4PhysicsTable* theLabTimeElectronTable ;
    static G4PhysicsTable* theLabTimePositronTable ;
    static G4PhysicsTable* theProperTimeElectronTable ;
    static G4PhysicsTable* theProperTimePositronTable ;

    //processes inherited from G4VeEnergyLoss 
    //register themselves  in the static array Recorder
    //for electrons/positrons separately
    //nb of contributing processes = NbOfProcesses
    static G4int NbOfProcesses;

    static G4int CounterOfElectronProcess;
    static G4int CounterOfPositronProcess ;          
    static G4PhysicsTable** RecorderOfElectronProcess;
    static G4PhysicsTable** RecorderOfPositronProcess;
    
  private:
     
    static G4int NbinEloss;               // number of bins in table,
                                          // calculated in BuildPhysicTable
    static G4double LowerBoundEloss;
    static G4double UpperBoundEloss;
    static G4double RTable,LOGRTable;    // LOGRTable=log(UpperBoundEloss-
                                         // LowerBoundEloss)/NbinEloss
                                         // RTable = exp(LOGRTable)

    //for interpolation within the tables
    static G4PhysicsTable* theeRangeCoeffATable;
    static G4PhysicsTable* theeRangeCoeffBTable;
    static G4PhysicsTable* theeRangeCoeffCTable;
    static G4PhysicsTable* thepRangeCoeffATable;
    static G4PhysicsTable* thepRangeCoeffBTable;
    static G4PhysicsTable* thepRangeCoeffCTable;
    

  public: // With description
     
};
 
#include "G4VeEnergyLoss.icc"

#endif
 
