// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VIMuEnergyLoss.hh,v 1.2 2000-06-07 16:50:21 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// -------------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4VIMuEnergyLoss physics process -----------
//                by Laszlo Urban, September 1997
// ********************************************************************
// It is the implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the continuous energy loss for muons.
// Processes giving contribution to the continuous loss :
//   ionisation (= cont.ion.loss + delta ray production)
//   bremsstrahlung
//   e+e- pair production
//   can be added more easily ..........
// This class creates static muplus/muminus dE/dx and range tables ,
// which tables can be used by other processes.
// ************************************************************
// some corrections by L.Urban on 27/05/98 , (but other corrections come soon!) 
// ------------------------------------------------------------
 
#ifndef G4VIMuEnergyLoss_h
#define G4VIMuEnergyLoss_h 1
 
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include "globals.hh"
#include "Randomize.hh"
#include "G4IVContinuousDiscreteProcess.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4EnergyLossTables.hh"
#include "G4VParticleChange.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
 
 
class G4VIMuEnergyLoss : public G4IVContinuousDiscreteProcess
 
{
  public:
    G4VIMuEnergyLoss(const G4String& );
    G4VIMuEnergyLoss(G4VIMuEnergyLoss &);

    virtual ~G4VIMuEnergyLoss();

    G4bool IsApplicable(const G4ParticleDefinition&);

  private:

  // hide  assignment operator 
    G4VIMuEnergyLoss & operator=(const G4VIMuEnergyLoss &right);

  public:

    G4double GetContinuousStepLimit(
                                    const G4Track& track,
                                    G4double previousStepSize,
                                    G4double currentMinimumStep,
                                    G4double& currentSafety) ; 


    G4VParticleChange* AlongStepDoIt(const G4Track& track ,const G4Step& Step) ;


    virtual G4VParticleChange* PostStepDoIt(const G4Track& track,const G4Step& Step) = 0 ;

  // Build energy loss table (total continuous energy loss)
  
    void BuildDEDXTable(const G4ParticleDefinition& aParticleType);

  //----------------------------------------------
  //  public functions .........................
  //  get the number of processes contributing to the cont.energy loss
    static G4int GetNUMBEROFPROCESSES()    { return NUMBEROFPROCESSES; }; 
  //  set the number of processes contributing to the cont.energy loss
    static void SetNUMBEROFPROCESSES(G4int number)
                                { NUMBEROFPROCESSES=number ; }; 
  //  Increment the number of processes contributing to the cont.energy loss
    static void PlusNUMBEROFPROCESSES()
                                { NUMBEROFPROCESSES++  ; }; 
  //  decrement the number of processes contributing to the cont.energy loss
    static void MinusNUMBEROFPROCESSES()
                                { NUMBEROFPROCESSES--  ; }; 


  //*****************************************************************************
  //

    G4double GetConstraints(const G4DynamicParticle *aParticle,
                            G4Material *aMaterial);
                                       
    G4double GetLossWithFluct(const G4DynamicParticle *aParticle,
                              G4Material *aMaterial) ;

  protected:

    G4PhysicsTable* theLossTable ;
   
    static G4PhysicsTable* theDEDXmuplusTable ;
    static G4PhysicsTable* theDEDXmuminusTable ;
    static G4PhysicsTable* theRangemuplusTable ;
    static G4PhysicsTable* theRangemuminusTable ;
    static G4PhysicsTable* theInverseRangemuplusTable ;
    static G4PhysicsTable* theInverseRangemuminusTable ;

    static G4PhysicsTable* theLabTimemuplusTable ;
    static G4PhysicsTable* theLabTimemuminusTable ;
    
    static G4PhysicsTable* theProperTimemuplusTable ;
    static G4PhysicsTable* theProperTimemuminusTable ;

    static G4PhysicsTable* themuplusRangeCoeffATable;
    static G4PhysicsTable* themuplusRangeCoeffBTable;
    static G4PhysicsTable* themuplusRangeCoeffCTable;
    static G4PhysicsTable* themuminusRangeCoeffATable;
    static G4PhysicsTable* themuminusRangeCoeffBTable;
    static G4PhysicsTable* themuminusRangeCoeffCTable;

    static G4double CutInmupluslossTable;
    static G4double CutInmuminuslossTable;

  //  processes inherited from G4VIMuEnergyLoss 
  //   register themselves  in the static array Recorder
  //  nb of contributing processes = NUMBEROFPROCESSES
    static G4int NUMBEROFPROCESSES ;  
    static G4PhysicsTable** RecorderOfmuplusProcess;
    static G4PhysicsTable** RecorderOfmuminusProcess;
    static G4int CounterOfmuplusProcess ;
    static G4int CounterOfmuminusProcess ;

    G4double RTable,LOGRTable; // LOGRTable=log(HighestKineticEnergy
                               //          /LowestKineticEnergy)/TotBin
                               //   RTable = exp(LOGRTable)

    // cut in range
    G4double CutInRange ;
    // last cut in range 
    G4double lastCutInRange ;

    // particle mass
    G4double ParticleMass;

    G4double BIGSTEP ;

  private:
  //  private functions ..................................

    void BuildRangeTable(const G4ParticleDefinition& aParticleType);

    void BuildInverseRangeTable(
                                const G4ParticleDefinition& aParticleType);

    void BuildTimeTables(const G4ParticleDefinition& aParticleType);

    void BuildLabTimeVector(G4int materialIndex,
                          G4PhysicsLogVector* rangeVector);

    void BuildProperTimeVector(G4int materialIndex,
                          G4PhysicsLogVector* rangeVector);

    void InvertRangeVector(G4int materialIndex,
                          G4PhysicsLogVector* rangeVector);


    void BuildRangeVector(G4int materialIndex,
                          G4PhysicsLogVector* rangeVector);

     G4double LabTimeIntLog(G4PhysicsVector* physicsVector,G4int nbin);
     G4double ProperTimeIntLog(G4PhysicsVector* physicsVector,G4int nbin);

    G4double RangeIntLin(G4PhysicsVector* physicsVector,G4int nbin);

    G4double RangeIntLog(G4PhysicsVector* physicsVector,G4int nbin);

    void BuildRangeCoeffATable(const G4ParticleDefinition& aParticleType);
    void BuildRangeCoeffBTable(const G4ParticleDefinition& aParticleType);
    void BuildRangeCoeffCTable(const G4ParticleDefinition& aParticleType);

  private:

  G4PhysicsTable* theDEDXTable;

  G4PhysicsTable* theRangeTable;
  G4PhysicsTable* theInverseRangeTable;

  G4PhysicsTable* theLabTimeTable;
  G4PhysicsTable* theProperTimeTable;

  G4PhysicsTable** RecorderOfProcess;
  G4int CounterOfProcess;
 
  //  private data members ...............................

    // fdEdx=(-dE/dx)
    // computed in GetConstraints at every call;
    G4double fdEdx;

    // fRangeNow is the actual range of the particle 
    //  computed in GetConstraints
    G4double fRangeNow ;

    // fMeanLoss is the energyloss without fluctuation
    //  computed in AlongStepDoIt ;
    G4double fMeanLoss ;

    // EnergyBinNumber,RangeCoeffA,... are needed to compute range
    G4int EnergyBinNumber ;
    G4double RangeCoeffA,RangeCoeffB,RangeCoeffC ;

   //................................................................
    G4PhysicsTable* theRangeCoeffATable;
    G4PhysicsTable* theRangeCoeffBTable;
    G4PhysicsTable* theRangeCoeffCTable;

    // dToverTini is the maximum allowed deltarange/range in one Step
    //  ( set in this class for the moment)
    const G4double dToverTini; 

    // LowestKineticEnergy = lower limit of particle kinetic energy
    // HighestKineticEnergy = upper limit of particle kinetic energy 
    // TotBin = number of bins 
    //  ---------in the energy loss/range tables-------------------
    const G4double LowestKineticEnergy;
    const G4double HighestKineticEnergy;
    G4int TotBin;// number of bins in table, calculated in BuildPhysicsTable
                 //  from LowestKineticEnergy,HighestKineticEnergy and  
                 //  dToverTini

    // variables for the integration routines
    G4double taulow,tauhigh,ltaulow,ltauhigh;

    //  cuts in kinetic energy ........ 
    G4double* ParticleCutInKineticEnergy ;
    G4double ParticleCutInKineticEnergyNow ; 

    // ...............
    const G4Electron* theElectron;
    const G4Positron* thePositron;
    const G4MuonPlus* theMuonPlus;
    const G4MuonMinus* theMuonMinus;

    // data members to speed up the fluctuation calculation
    G4Material *lastMaterial ;
    G4double f1Fluct,f2Fluct,e1Fluct,e2Fluct,rateFluct,ipotFluct;
    G4double e1LogFluct,e2LogFluct,ipotLogFluct;
    const G4double MaxExcitationNumber ;
    const G4double probLimFluct ;
    const long nmaxDirectFluct,nmaxCont1,nmaxCont2 ;

};
 
#include "G4VIMuEnergyLoss.icc"

#endif
