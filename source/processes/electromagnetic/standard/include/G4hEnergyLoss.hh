// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hEnergyLoss.hh,v 1.8 2000-02-22 10:37:49 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hEnergyLoss physics process -----------
//                by Laszlo Urban, 30 May 1997 
//
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
//
// **************************************************************************** 
// It is assumed that the cut in range is the same for all the charged hadrons! 
// ****************************************************************************
//
// 7/10/98 some bugs fixed + some cleanup , L.Urban 
// 22/10/98 cleanup , L.Urban
// 02/02/99 several bugs fixed, L.Urban
// 10/02/00  modifications , new e.m. structure, L.Urban
//

#ifndef G4hEnergyLoss_h
#define G4hEnergyLoss_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VEnergyLoss.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Electron.hh"
#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
 
class G4EnergyLossMessenger;
 
class G4hEnergyLoss : public G4VEnergyLoss
 
{
  public:

    G4hEnergyLoss(const G4String& );

    ~G4hEnergyLoss();

    G4bool IsApplicable(const G4ParticleDefinition&);

    G4double GetContinuousStepLimit(
                                    const G4Track& track,
                                    G4double previousStepSize,
                                    G4double currentMinimumStep,
                                    G4double& currentSafety) ; 

    G4VParticleChange* AlongStepDoIt(const G4Track& track ,const G4Step& Step) ;

    virtual G4double GetMeanFreePath(
                                      const G4Track& track,
                                      G4double previousStepSize,
                                      G4ForceCondition* condition
                                                         ) = 0 ;

    virtual G4VParticleChange* PostStepDoIt(const G4Track& track,
                                               const G4Step& Step) = 0 ;


  protected:

    virtual G4double GetConstraints(const G4DynamicParticle *aParticle,
                            G4Material *aMaterial);
                                       
  private:

  // hide  assignment operator 

    G4hEnergyLoss(G4hEnergyLoss &);
    G4hEnergyLoss & operator=(const G4hEnergyLoss &right);


// =====================================================================

  public:


// ====================================================================
//  static part of the class

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
    static void SetNbinEloss(G4int nb)           {NbinEloss=nb;};

    static G4double GetLowerBoundEloss() {return LowerBoundEloss;};
    static G4double GetUpperBoundEloss() {return UpperBoundEloss;};
    static G4int    GetNbinEloss()       {return NbinEloss;};

  protected:

    static void BuildDEDXTable(const G4ParticleDefinition& aParticleType);

// ====================================================================

  public:

  protected:

    static G4PhysicsTable* theDEDXpTable ;
    static G4PhysicsTable* theDEDXpbarTable ;
    static G4PhysicsTable* theRangepTable ;
    static G4PhysicsTable* theRangepbarTable ;

    //inverse of the range tables
    static G4PhysicsTable* theInverseRangepTable ;
    static G4PhysicsTable* theInverseRangepbarTable ;

    //lab and proper time tables
    static G4PhysicsTable* theLabTimepTable ;
    static G4PhysicsTable* theLabTimepbarTable ;

    static G4PhysicsTable* theProperTimepTable ;
    static G4PhysicsTable* theProperTimepbarTable ;

    //  processes inherited from G4hEnergyLoss 
    //   register themselves  in the static array Recorder

    static G4int NbOfProcesses     ;
    static G4PhysicsTable** RecorderOfpProcess;
    static G4PhysicsTable** RecorderOfpbarProcess;
    static G4int CounterOfpProcess ;
    static G4int CounterOfpbarProcess ;

    // cut in range
    static G4double ptableElectronCutInRange;
    static G4double pbartableElectronCutInRange;

    static G4double Charge ;

    G4PhysicsTable* theLossTable ;

    G4double linLossLimit ;
   
    G4double MinKineticEnergy ;

    G4double fdEdx;      // computed in GetContraints
    G4double fRangeNow ; // computed in GetContraints

    static const G4Proton* theProton ;
    static const G4AntiProton* theAntiProton ;

  private:

    static G4int NbinEloss;               // number of bins in table,
                                          // calculated in BuildPhysicTable
    static G4double LowerBoundEloss;
    static G4double UpperBoundEloss;
    static G4double RTable,LOGRTable;    // LOGRTable=log(UpperBoundEloss-
                                         // LowerBoundEloss)/NbinEloss
                                         // RTable = exp(LOGRTable)

    static G4PhysicsTable** RecorderOfProcess;
    static G4int CounterOfProcess;


    static G4PhysicsTable* thepRangeCoeffATable;
    static G4PhysicsTable* thepRangeCoeffBTable;
    static G4PhysicsTable* thepRangeCoeffCTable;
    static G4PhysicsTable* thepbarRangeCoeffATable;
    static G4PhysicsTable* thepbarRangeCoeffBTable;
    static G4PhysicsTable* thepbarRangeCoeffCTable;

    static G4PhysicsTable* theDEDXTable ;


};
 
#include "G4hEnergyLoss.icc"

#endif
 


