// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IMultipleScattering.hh,v 1.2 1999-12-15 14:51:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id:
// --------------------------------------------------------------
//    GEANT 4 class header file
//  
//    For information related to this code contact:
//    CERN, IT Division, ASD Group
//    History: based on object model of
//    2nd December 1995, G.Cosmo
//    --------- G4IMultipleScattering physics process --------
//               by Laszlo Urban, October 1997
// **************************************************************
//      UNIVERSAL: for arbitrary single charged particle
//  09/12/98:      charge can be != +- 1 !!!!   L.Urban
// --------------------------------------------------------------
// *****************************************************************
// It is the first implementation of the multiple scattering process
//   using an INTEGRAL APPROACH instead of the differential
//   one used in the standard implementation .
// *****************************************************************
//                by Laszlo Urban, 23 June 1998
// -----------------------------------------------------------------
// 27/10/98: cleanup , L.Urban

#ifndef G4IMultipleScattering_h
#define G4IMultipleScattering_h 1

#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "globals.hh"
#include "Randomize.hh"
#include "G4EnergyLossTables.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4Proton.hh"
#include "G4PhysicsLogVector.hh"
#include "G4GPILSelection.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4ElementTable.hh"
#include "G4ElementVector.hh"
#include "G4VParticleChange.hh"

class G4IMultipleScattering : public G4VContinuousDiscreteProcess

{
 public:

   G4IMultipleScattering(const G4String& processName="Imsc") ;

   ~G4IMultipleScattering() ;
          
   G4bool IsApplicable ( const G4ParticleDefinition& ) ;

   void SetPhysicsTableBining(G4double lowE,G4double highE,G4int nBins);
 
   void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) ;

   void BuildIntegralITable(const G4ParticleDefinition& aParticleType) ;
   void BuildIntegralJTable(const G4ParticleDefinition& aParticleType) ;

   G4double GetIntegralI(const G4ParticleDefinition *aParticle,
                          G4double KineticEnergy,G4Material* aMaterial) ;
   G4double GetIntegralJ(const G4ParticleDefinition *aParticle,
                           G4double KineticEnergy,G4Material* aMaterial) ;

   void PrintInfoDefinition();

   G4double GetContinuousStepLimit(const G4Track& aTrack,
                                   G4double previousStepSize,
                                   G4double currentMinimumStep,
                                   G4double& currentSafety) ; 

   G4double GetMeanFreePath(const G4Track& aTrack,
                            G4double previousStepSize,
                            G4ForceCondition* condition) ;

   G4VParticleChange* AlongStepDoIt(const G4Track& aTrack,const G4Step& aStep) ;

   G4VParticleChange* PostStepDoIt(const G4Track& aTrack,const G4Step& aStep) ;

 protected:

   G4double ComputeTransportCrossSection(
                             const G4ParticleDefinition& aParticleType,
                                   G4double KineticEnergy,
                                   G4double AtomicNumber) ;

   G4double TrueToGeomTransformation(const G4DynamicParticle* aParticle,
                                     G4Material* aMaterial,
                                     G4double truePathLength) ;

   
 private:

 //  hide assignment operator as  private
   G4IMultipleScattering & operator = (const G4IMultipleScattering &right) ;
   G4IMultipleScattering ( const G4IMultipleScattering &) ;


 private:

  // data members ...................................................
   G4PhysicsTable* theTransportMeanFreePathTable;
   G4PhysicsTable* theIntegralITable ;
   G4PhysicsTable* theIntegralJTable ;

   G4double fTransportMeanFreePath ;

   G4double fMeanLateralDisplacement ;

   G4double LowestKineticEnergy ;
   G4double HighestKineticEnergy ;
   G4int TotBin ;

   const G4Electron* theElectron ;
   const G4Positron* thePositron ;
   
   G4Material* lastMaterial;
   G4double lastKineticEnergy;

   // GeomStepFinal is the geom.steplength at the end of the AlongStep loop
   // the others are some 'cache' variables 
   G4double GeomStepFinal ;
   G4double tLast,zLast,CosTheta ;

   G4double biglambda ;

   //parameters for low energy extrapolation of dE/dx and lambda
   const G4double plowloss,plowlambda ;
  
   G4int NumberOfBuildPhysicsTableCalls ;

   G4double tuning ;
};

#include "G4IMultipleScattering.icc"

#endif
 

 
