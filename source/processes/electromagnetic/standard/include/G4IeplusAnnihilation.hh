// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IeplusAnnihilation.hh,v 1.1 1999-01-07 16:11:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id: 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4IeplusAnnihilation process ------
//                   by Michel Maire, 7 july 1996
// ************************************************************
// ************************************************************
// It is the first implementation of the
//          eplusANNIHILATION  PROCESS
//   using an INTEGRAL APPROACH instead of the differential
//   one used in the standard implementation .
// ************************************************************
//                by Laszlo Urban, 23 June 1998
// ----------------------------------------------------------
// 28/10/98: cleanup   L. Urban


#ifndef G4IeplusAnnihilation_h
#define G4IeplusAnnihilation_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4IVRestDiscreteProcess.hh"
#include "G4EnergyLossTables.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh" 
#include "G4PhysicsLinearVector.hh" 
#include "G4ElementTable.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Step.hh"

class G4IeplusAnnihilation : public G4IVRestDiscreteProcess
 
{ 
  public:
 
     G4IeplusAnnihilation(const G4String& processName ="Iannihil");
 
    ~G4IeplusAnnihilation();

     G4bool IsApplicable(const G4ParticleDefinition&);

     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);

     void PrintInfoDefinition();

     void BuildPhysicsTable(const G4ParticleDefinition& PositronType);

     G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                              G4double previousStepSize,
                                              G4ForceCondition* condition) ;   

     G4double GetMicroscopicCrossSection(G4DynamicParticle* aDynamicPositron,
                                         G4Element*         anElement);

     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step& aStep); 

     G4double GetMeanLifeTime(const G4Track& aTrack,
                              G4ForceCondition* condition);

     G4VParticleChange* AtRestDoIt(const G4Track& aTrack,
                                  const G4Step& aStep); 

  protected:

     virtual G4double ComputeMicroscopicCrossSection(G4double PositKinEnergy,
                                                     G4double AtomicNumber);

     virtual G4double ComputeMeanFreePath(G4double PositKinEnergy, 
                                          G4Material* aMaterial);


  private:
               
  // hide assignment operator as private 
      G4IeplusAnnihilation& operator=(const G4IeplusAnnihilation &right);
      G4IeplusAnnihilation(const G4IeplusAnnihilation& );

    void BuildNlambdaTable(const G4ParticleDefinition& aParticleType) ;
    void BuildNlambdaVector(const G4ParticleDefinition& aParticleType,
                                  G4int materialIndex,
                                  G4PhysicsLogVector* nlambdaVector) ;
    void BuildInverseNlambdaTable(
                           const G4ParticleDefinition& aParticleType) ;
    void InvertNlambdaVector(const G4ParticleDefinition& aParticleType,
                                   G4int materialIndex,
                                   G4PhysicsLogVector* nlambdaVector) ;

    void BuildCoeffATable(const G4ParticleDefinition& aParticleType) ;
    void BuildCoeffBTable(const G4ParticleDefinition& aParticleType) ;
    void BuildCoeffCTable(const G4ParticleDefinition& aParticleType) ;

    void TestOfInversion(const G4ParticleDefinition& aParticleType,
                                     G4int printflag) ;

  private:

    G4PhysicsTable*   theCrossSectionTable;    // table for crossection

    G4PhysicsTable* theMeanFreePathTable ;

    G4PhysicsTable* theNlambdaTable;
    G4PhysicsTable* theInverseNlambdaTable;

    G4PhysicsTable* theCoeffATable;
    G4PhysicsTable* theCoeffBTable;
    G4PhysicsTable* theCoeffCTable;

    G4double LowestEnergyLimit ;      // low  energy limit of the crossection formula
    G4double HighestEnergyLimit ;     // high energy limit of the crossection formula 
    G4int NumbBinTable ;              // number of bins in the crossection table
    G4int NumberOfBuildPhysicsTableCalls ;    
     
    G4double LowestKineticEnergy ;
    G4double HighestKineticEnergy;
    G4int TotBin ;

    G4double RTable ;

};

#include "G4IeplusAnnihilation.icc"
  
#endif
 
