// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAIonisation.hh,v 1.7 2000-09-22 14:44:07 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4hIonisation physics process -----------
//                by V. Grichine, 30 Nov 1997 
// ************************************************************
// It is the first implementation of the NEW IONISATION     
// PROCESS. ( delta rays + continuous energy loss)
// It calculates the ionisation for charged hadrons.      
// ************************************************************
//
// corrected by V. Grichine on 24/11/97
// 12.07.00 V.Grichine GetFreePath and GetdEdx were added
// 22.09.00 V.Grichine and K.Assamagan new function GetEnergyTransfer(G4d,G4i)

 
#ifndef G4PAIonisation_h
#define G4PAIonisation_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VPAIenergyLoss.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Electron.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsVector.hh"
 
 
class G4PAIonisation : public G4VPAIenergyLoss 
 
{
  public:
 
     G4PAIonisation( const G4String& materialName,
                     const G4String& processName = "PAIonisation"); 

    ~G4PAIonisation();

     G4bool IsApplicable(const G4ParticleDefinition&);

  private:

  // hide assignment operator 

    G4PAIonisation & operator=(const G4PAIonisation &right);
    G4PAIonisation(const G4PAIonisation&);

  public:

  // post Step functions ....................................... 


    G4double GetConstraints(const G4DynamicParticle *aParticle,
                            G4Material *aMaterial);
 
     G4VParticleChange *PostStepDoIt( const G4Track& track,         
                                      const G4Step& Step      ) ;                 


     //  void BuildLossTable(const G4ParticleDefinition& aParticleType);

     void BuildLambdaTable(const G4ParticleDefinition& aParticleType);

     void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);

     void  BuildPAIonisationTable() ;


     virtual G4double ComputeMicroscopicCrossSection(
                            const G4ParticleDefinition& aParticleType,
                            G4double KineticEnergy,
                            G4double AtomicNumber);

  // Along step DoIt stuff

    G4double GetContinuousStepLimit(
                                    const G4Track& track,
                                    G4double previousStepSize,
                                    G4double currentMinimumStep,
                                    G4double& currentSafety) ; 

    G4double GetMeanFreePath( const G4Track& track,
                               G4double previousStepSize,
                               G4ForceCondition* condition   ) ;

    G4double GetFreePath( G4double scaledTkin, G4double charge2 ) ;

    G4double GetdEdx( G4double scaledTkin, G4double charge2 ) ;

    G4VParticleChange* AlongStepDoIt(const G4Track& track ,const G4Step& Step) ;

    G4double GetLossWithFluct(G4double Step,
                              const G4DynamicParticle *aParticle,
                              G4Material *aMaterial) ;

    G4double GetRandomEnergyTransfer( G4double scaledTkin ) ;
    G4double GetEnergyTransfer( G4int iPlace, G4double position, G4int iTransfer  ) ;

  //    static                                   

    static G4double GetMaxKineticEnergy() { return HighestKineticEnergy ; } ;
    static G4double GetMinKineticEnergy() { return  LowestKineticEnergy ; } ;
    static G4int    GetBinNumber()        { return               TotBin ; } ;

  // Access to Sandia table coefficients

    G4int    GetSandiaIntervalNumber()    const {return fSandiaIntervalNumber;};
    G4double GetSandiaPhotoAbsCof(G4int, G4int) const;

    // Compute Sandia photoabsorption coefficient matrix
    void ComputeSandiaPhotoAbsCof() ;
 

  private:

    //  private data members 

    G4int fMatIndex ;  // index of material, where dE/dx is calculated

    G4PhysicsTable* theMeanFreePathTable;

    // LowestKineticEnergy = lower limit of particle kinetic energy
    // HighestKineticEnergy = upper limit of particle kinetic energy 
    // TotBin = number of bins in the energy ionisation loss table 
    //

    static const G4double LowestKineticEnergy;
    static const G4double HighestKineticEnergy;
    static G4int TotBin;
    static G4PhysicsLogVector* fProtonEnergyVector ;

    // cut in range

    G4double CutInRange ;
    G4double lastCutInRange ;

    // particles , cuts in kinetic energy 

    const G4Electron* theElectron;
    const G4Proton* theProton;
    const G4AntiProton* theAntiProton;

    const G4double* ParticleCutInKineticEnergy;
    const G4double* DeltaCutInKineticEnergy ; 
 
    G4double ParticleCutInKineticEnergyNow ; 
    G4double DeltaCutInKineticEnergyNow ;

    // Sandia Photo Absorption Coefficients
       
    G4double** fSandiaPhotoAbsCof ; 
    G4int      fSandiaIntervalNumber ;


};
 
#include "G4PAIonisation.icc"
 
#endif
 







