// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GammaConversion.hh,v 1.2 1999-12-15 14:51:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4GammaConversion physics process ------
//                   by Michel Maire, 24 May 1996
// ************************************************************
// 11-06-96, Added GetRandomAtom() method and new data member
//           for cumulative total cross section, by M.Maire
// 21-06-96, SetCuts inplementation, M.Maire
// 16-09-96, Dynamical array PartialSumSigma, M.Maire
// 14-01-97, crossection table + meanfreepath table.
//           PartialSumSigma removed, M.Maire
// 14-03-97, new physics scheme for geant4alpha, M.Maire
// 13-08-98, new methods SetBining() PrintInfo() 
// ------------------------------------------------------------

#ifndef G4GammaConversion_h
#define G4GammaConversion_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Element.hh"
#include "G4Gamma.hh" 
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Step.hh"
 
class G4GammaConversion : public G4VDiscreteProcess
 
{  
  public:
 
     G4GammaConversion(const G4String& processName ="conv");
 
    ~G4GammaConversion();

     G4bool IsApplicable(const G4ParticleDefinition&);
     
     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);

     void BuildPhysicsTable(const G4ParticleDefinition& GammaType);
     
     void PrintInfoDefinition();

     G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double previousStepSize,
                              G4ForceCondition* condition);

     G4double GetMicroscopicCrossSection(const G4DynamicParticle* aDynamicGamma,
                                         G4Element*         anElement);

     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step& aStep);
 
  protected:

     virtual G4double ComputeMicroscopicCrossSection(G4double GammaEnergy, 
                                             G4double AtomicNumber);

     virtual G4double ComputeMeanFreePath (G4double GammaEnergy, 
                                           G4Material* aMaterial);

  private:

     G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicGamma,
                                 G4Material* aMaterial);

     G4double ScreenFunction1(G4double ScreenVariable);

     G4double ScreenFunction2(G4double ScreenVariable);
     
  private:
  
     // hide assignment operator as private 
     G4GammaConversion& operator=(const G4GammaConversion &right);
     G4GammaConversion(const G4GammaConversion& );
     
  private:
  
     G4PhysicsTable* theCrossSectionTable;    // table for crossection
     G4PhysicsTable* theMeanFreePathTable;
     
     G4double LowestEnergyLimit ;      // low  energy limit of the crossection formula
     G4double HighestEnergyLimit ;     // high energy limit of the crossection formula 
     G4int NumbBinTable ;              // number of bins in the crossection table

     G4double MeanFreePath;            // actual Mean Free Path (current medium)
};

#include "G4GammaConversion.icc"
  
#endif
 
