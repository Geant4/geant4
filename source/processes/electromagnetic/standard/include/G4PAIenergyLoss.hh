// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAIenergyLoss.hh,v 1.6 2000-02-22 10:37:49 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4PAIenergyLoss physics process -----------
//                by V. Grichine, 30 Nov 97 
// ************************************************************             
// It is the first implementation of the NEW UNIFIED ENERGY LOSS PROCESS.
// It calculates the continuous energy loss for charged hadrons.
// Processes giving contribution to the continuous loss :
//   ionisation (= cont.ion.loss + delta ray production)
//   can be added more easily ..........
// This class creates static proton/antiproton dE/dx and range tables ,
// which tables can be used by other processes.
// The energy loss for other charged hadrons is calculated from the p/pbar
// tables with scaled kinetic energy.
//******************************************************************************* 
//* It is assumed that the cut in range is the same for all the charged hadrons!* 
//*******************************************************************************
// corrected by V. Grichine on 24/11/97 
// corrected by L. Urban    on 27/05/98  (other corrections come soon!)
// 10/02/00  modifications , new e.m. structure, L.Urban
// ------------------------------------------------------------
 
#ifndef G4PAIenergyLoss_h
#define G4PAIenergyLoss_h 1
 
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include "g4rw/tpordvec.h"

#include "globals.hh"
#include "Randomize.hh"
#include "G4VEnergyLoss.hh"
#include "G4MaterialTable.hh"
#include "G4ElementTable.hh"
#include "G4ElementVector.hh"
#include "G4VParticleChange.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsFreeVector.hh"
 
 
class G4PAIenergyLoss : public G4VEnergyLoss
 
{
  public:
    G4PAIenergyLoss(const G4String& );
    G4PAIenergyLoss(G4PAIenergyLoss &);

    ~G4PAIenergyLoss();

    G4bool IsApplicable(const G4ParticleDefinition&);

  private:

  // hide  assignment operator 
    G4PAIenergyLoss & operator=(const G4PAIenergyLoss &right);

  public:

    virtual G4double 
    GetContinuousStepLimit( const G4Track& track,
                                  G4double previousStepSize,
                                  G4double currentMinimumStep,
                                  G4double& currentSafety      )= 0 ; 


     /* *************************
    G4VParticleChange* AlongStepDoIt(const G4Track& track ,const G4Step& Step)
{ 
//  clear  NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();
    return pParticleChange;
} ;
      **************************** */

    virtual G4double GetMeanFreePath(
                                      const G4Track& track,
                                      G4double previousStepSize,
                                      G4ForceCondition* condition
                                                         ) = 0 ;

    virtual G4VParticleChange* PostStepDoIt(const G4Track& track,const G4Step& Step) = 0 ;

  // Build energy loss table (total continuous energy loss)
  
    static void BuildDEDXTable(const G4ParticleDefinition& aParticleType);

  //----------------------------------------------
  //  public functions .........................
  //
     G4double GetcurrentInteractionLength() const { return currentInteractionLength; } ;

  // returns the actual dE/dx in internal GEANT4 units
    G4double GetdEdx() const { return fdEdx; };

  // returns the actual range of the particle
    G4double GetRangeNow() const { return fRangeNow; } ;

  // returns the mean energy loss (i.e. the loss
  //  without fluctuations
    G4double GetMeanLoss() const { return fMeanLoss; } ;  


    G4double GetConstraints(const G4DynamicParticle *aParticle,
                            G4Material *aMaterial);

  // static 
    G4PhysicsTable* GetPAItransferBank(){ return fPAItransferBank ; } ;

    static G4double GetMaxKineticEnergy() { return UpperBoundEloss ; } ;
    static G4double GetMinKineticEnergy() { return  LowerBoundEloss ; } ;
    static G4int    GetBinNumber()        { return               NbinEloss ; } ; 

  public:  // With description

    static void  SetNbOfProcesses(G4int nb) {NbOfProcesses=nb;};
    // Sets number of processes giving contribution to the energy loss

    static void  PlusNbOfProcesses()        {NbOfProcesses++ ;};
    // Increases number of processes giving contribution to the energy loss

    static void  MinusNbOfProcesses()       {NbOfProcesses-- ;};
    // Decreases number of processes giving contribution to the energy loss

    static G4int GetNbOfProcesses()         {return NbOfProcesses;};
    // Gets number of processes giving contribution to the energy loss
    // ( default value = 1)

    static void SetLowerBoundEloss(G4double val) {LowerBoundEloss=val;};
    static void SetUpperBoundEloss(G4double val) {UpperBoundEloss=val;};
    static void SetNbinEloss(G4int nb)		 {NbinEloss=nb;};

    static G4double GetLowerBoundEloss() {return LowerBoundEloss;};
    static G4double GetUpperBoundEloss() {return UpperBoundEloss;};
    static G4int    GetNbinEloss()	 {return NbinEloss;};


  protected:

    // fRangeNow is the actual range of the particle 
    //  computed in GetConstraints
    G4double fRangeNow ;

    // fMeanLoss is the energyloss without fluctuation
    //  computed in AlongStepDoIt ;
    G4double fMeanLoss ;

    G4PhysicsTable* theLossTable ;
   
    static G4PhysicsTable* theDEDXpTable ;
    static G4PhysicsTable* theDEDXpbarTable ;
    static G4PhysicsTable* theRangepTable ;
    static G4PhysicsTable* theRangepbarTable ;
    static G4PhysicsTable* theInverseRangepTable ;
    static G4PhysicsTable* theInverseRangepbarTable ;
    static G4PhysicsTable* theLabTimepTable ;
    static G4PhysicsTable* theLabTimepbarTable ;
    static G4PhysicsTable* theProperTimepTable ;
    static G4PhysicsTable* theProperTimepbarTable ;

    //  bank of PAI energy transfer data

    // static 
    G4PhysicsTable* fPAItransferBank ;



  //  processes inherited from G4hEnergyLoss 
  //   register themselves  in the static array Recorder

    static G4int NbOfProcesses     ;
    static G4PhysicsTable** RecorderOfpProcess;
    static G4PhysicsTable** RecorderOfpbarProcess;
    static G4int CounterOfpProcess ;
    static G4int CounterOfpbarProcess ;

  private:
  //  private functions ..................................

  private:

  //  private data members ...............................

    static G4PhysicsTable* theDEDXTable ; 

    // fdEdx=(-dE/dx)
    // computed in GetConstraints at every call;
    // it can be used by other processes ( Cherenkov, ...)
    G4double fdEdx;

    G4double dToverTini ;


    // EnergyBinNumber,RangeCoeffA,... are needed to compute range
    G4int EnergyBinNumber ;
    G4double RangeCoeffA,RangeCoeffB,RangeCoeffC ;

   //................................................................
    static G4PhysicsTable* thepRangeCoeffATable;
    static G4PhysicsTable* thepRangeCoeffBTable;
    static G4PhysicsTable* thepRangeCoeffCTable;
    static G4PhysicsTable* thepbarRangeCoeffATable;
    static G4PhysicsTable* thepbarRangeCoeffBTable;
    static G4PhysicsTable* thepbarRangeCoeffCTable;


    // LowerBoundEloss = lower limit of particle kinetic energy
    // UpperBoundEloss = upper limit of particle kinetic energy 
    // NbinEloss = number of bins 
    //  ---------in the energy loss/range tables-------------------

    static G4double LowerBoundEloss;
    static G4double UpperBoundEloss;    
    static G4int NbinEloss;

    static G4double RTable,LOGRTable; // LOGRTable=log(UpperBoundEloss
                               //          /LowerBoundEloss)/NbinEloss
                               //   RTable = exp(LOGRTable)


    // cut in range
    static G4double CutInRange;

    //  cuts in kinetic energy ........ 
    G4double* ParticleCutInKineticEnergy ;
    G4double ParticleCutInKineticEnergyNow ; 

    // ...............
    const G4Electron* theElectron;
    const G4Proton* theProton;
    const G4AntiProton* theAntiProton;


};
 
#include "G4PAIenergyLoss.icc"

#endif
 









