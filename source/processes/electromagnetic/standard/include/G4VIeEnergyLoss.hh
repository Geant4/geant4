// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VIeEnergyLoss.hh,v 1.2 2000-06-07 17:00:48 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4VIeEnergyLoss physics process -----------
//                by Laszlo Urban, 20 March 1997 
// ************************************************************             
// It is the first implementation of the new unified Energy Loss process.
// It calculates the continuous energy loss for e+/e-.
// Processes giving contribution to the continuous loss :
//   ionisation (= cont.ion.loss + delta ray production)
//   bremsstrahlung (= cont.loss due to sooft brems+discrete bremsstrahlung)
//   can be added more easily ..........
// This class creates static dE/dx and range tables for e+ and e-,
// which tables can be used by other processes.
// ------------------------------------------------------------
//
//  27.05.98 OldGetRange removed + other corrs , L.Urban
//  26.10.98 revision , L.Urban
// ------------------------------------------------------------
 
#ifndef G4VIeEnergyLoss_h
#define G4VIeEnergyLoss_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4IVContinuousDiscreteProcess.hh"
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

class G4EnergyLossMessenger;
 
 
class G4VIeEnergyLoss : public G4IVContinuousDiscreteProcess
 
{
  public:
 
    G4VIeEnergyLoss(const G4String& );

    virtual ~G4VIeEnergyLoss();

    G4bool IsApplicable(const G4ParticleDefinition&);


  public:
  
    void BuildDEDXTable(const G4ParticleDefinition& aParticleType);

    G4double GetContinuousStepLimit(const G4Track& track,
                                    G4double previousStepSize,
                                    G4double currentMinimumStep,
                                    G4double& currentSafety);

    G4VParticleChange* AlongStepDoIt(const G4Track& track,
                                     const G4Step& Step) ;

    virtual G4VParticleChange* PostStepDoIt(const G4Track& track,
                                            const G4Step& Step) = 0;
                                            
                                            
  private:

    void BuildRangeTable(const G4ParticleDefinition& aParticleType);

    void BuildInverseRangeTable(const G4ParticleDefinition& aParticleType);

    void BuildTimeTables(const G4ParticleDefinition& aParticleType);

    void BuildRangeVector(G4int materialIndex,
                          G4PhysicsLogVector* rangeVector);

    void BuildLabTimeVector(G4int materialIndex,
                          G4PhysicsLogVector* rangeVector);

    void BuildProperTimeVector(G4int materialIndex,
                          G4PhysicsLogVector* rangeVector);

    void InvertRangeVector(G4int materialIndex,
                          G4PhysicsLogVector* rangeVector);

    G4double RangeIntLin(G4PhysicsVector* physicsVector,G4int nbin);

    G4double RangeIntLog(G4PhysicsVector* physicsVector,G4int nbin);
    G4double LabTimeIntLog(G4PhysicsVector* physicsVector,G4int nbin);
    G4double ProperTimeIntLog(G4PhysicsVector* physicsVector,G4int nbin);

    void BuildRangeCoeffATable(const G4ParticleDefinition& aParticleType);
    void BuildRangeCoeffBTable(const G4ParticleDefinition& aParticleType);
    void BuildRangeCoeffCTable(const G4ParticleDefinition& aParticleType);
    
    G4double GetConstraints(const G4DynamicParticle* aParticle,
                            G4Material* aMaterial); 
                                                                  
    G4double GetLossWithFluct(const G4DynamicParticle* aParticle,
                              G4Material* aMaterial,
                              G4double   threshold);

    // hide  assignment operator
    G4VIeEnergyLoss (G4VIeEnergyLoss &); 
    G4VIeEnergyLoss & operator=(const G4VIeEnergyLoss &right);

  protected:

    G4PhysicsTable* theLossTable;
     
    G4double ParticleMass;               // heavily used

  private:

    G4PhysicsTable* theDEDXTable;
 
    G4PhysicsTable* theRangeTable;
    G4PhysicsTable* theInverseRangeTable;

    G4PhysicsTable* theLabTimeTable ;
    G4PhysicsTable* theProperTimeTable ;
  
    G4int            CounterOfProcess;
    G4PhysicsTable** RecorderOfProcess;
                                            
    G4double fdEdx;                       // computed in GetConstraints
    G4double fRangeNow;                   // computed in GetConstraints

    G4int    EnergyBinNumber;             // computed in GetConstraints
                                          // (needed to compute range)   
    G4int TotBin;                         // number of bins in table, 
                                          // calculated in BuildPhysicTable
    G4double LowestKineticEnergy;
    G4double HighestKineticEnergy;
    G4double RTable,LOGRTable;           // LOGRTable=log(HighestKineticEnergy-
                                         // LowestKineticEnergy)/TotBin
                                         // RTable = exp(LOGRTable)

    G4PhysicsTable* theRangeCoeffATable;
    G4PhysicsTable* theRangeCoeffBTable;
    G4PhysicsTable* theRangeCoeffCTable;
    
    // variables for the integration routines
    G4double taulow,tauhigh,ltaulow,ltauhigh;
    
    // data members to speed up the fluctuation calculation
    G4Material* lastMaterial;
    G4int imat;
    G4double f1Fluct,f2Fluct,e1Fluct,e2Fluct,rateFluct,ipotFluct;
    G4double e1LogFluct,e2LogFluct,ipotLogFluct;
    const G4double MaxExcitationNumber ;
    const G4double probLimFluct ;
    const long nmaxDirectFluct,nmaxCont1,nmaxCont2 ;
    
 //  
 // static part of the class
 //
 
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

    //processes inherited from G4VIeEnergyLoss 
    //register themselves  in the static array Recorder
    //for electrons/positrons separately
    //nb of contributing processes = NbOfProcesses
    static G4int NbOfProcesses;
    static G4int CounterOfElectronProcess;
    static G4int CounterOfPositronProcess ;          
    static G4PhysicsTable** RecorderOfElectronProcess;
    static G4PhysicsTable** RecorderOfPositronProcess;
    
  private:
     
    //for interpolation within the tables
    static G4PhysicsTable* theeRangeCoeffATable;
    static G4PhysicsTable* theeRangeCoeffBTable;
    static G4PhysicsTable* theeRangeCoeffCTable;
    static G4PhysicsTable* thepRangeCoeffATable;
    static G4PhysicsTable* thepRangeCoeffBTable;
    static G4PhysicsTable* thepRangeCoeffCTable;
    
    static G4double dRoverRange;     // dRoverRange is the maximum allowed
                                     // deltarange/range in one Step 
    static G4double finalRange;      // final step before stopping
    
    static G4bool   rndmStepFlag;    // control the randomization of the step
    static G4bool   EnlossFlucFlag;  // control the energy loss fluctuation
    
    static G4EnergyLossMessenger* eLossMessenger;
         
  public:
     
    static void  SetNbOfProcesses(G4int nb) {NbOfProcesses=nb;};
    static void  PlusNbOfProcesses()        {NbOfProcesses++ ;};
    static void  MinusNbOfProcesses()       {NbOfProcesses-- ;};                                      
    static G4int GetNbOfProcesses()         {return NbOfProcesses;};
    
    static void SetRndmStep     (G4bool   value) {rndmStepFlag   = value;}
    static void SetEnlossFluc   (G4bool   value) {EnlossFlucFlag = value;}
    static void SetStepFunction (G4double c1, G4double c2) 
                                      {dRoverRange = c1; finalRange = c2;}         
};
 
#include "G4VIeEnergyLoss.icc"

#endif
 
