// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ComptonScattering.hh,v 1.4 1999-12-17 16:20:49 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ------------ G4ComptonScattering physics process ------
//                   by Michel Maire, April 1996
//
// 10-06-96, updated by M.Maire 
// 21-06-96, SetCuts implementation, M.Maire
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 17-02-97, New Physics scheme
// 25-02-97, GetMeanFreePath() now is public function
// 12-03-97, new physics scheme again
// 13-08-98, new methods SetBining()  PrintInfo() 
// ------------------------------------------------------------

// class description
//
// inherit from G4VDiscreteProcess
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4ComptonScattering_h
#define G4ComptonScattering_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh" 
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ComptonScattering : public G4VDiscreteProcess
 
{ 
  public:  // with description
 
     G4ComptonScattering(const G4String& processName ="compt");
 
    ~G4ComptonScattering();

     G4bool IsApplicable(const G4ParticleDefinition&);
       // true for Gamma only.
     
     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);
       // Allows to define the binning of the PhysicsTables, 
       // before to build them.
     
     void BuildPhysicsTable(const G4ParticleDefinition& GammaType);
       // It builds the total CrossSectionPerAtom table, for Gamma,
       // and for every element contained in the elementTable.
       // It builds the MeanFreePath table, for Gamma,
       // and for every material contained in the materialTable.       
       // This function overloads a virtual function of the base class.
       // It is invoked by the G4ParticleWithCuts::SetCut() method. 
     
     void PrintInfoDefinition();
       // Print few lines of informations about the process: validity range,
       // origine ..etc..
       // Invoked by BuildPhysicsTable(). 
     
     G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double previousStepSize,
                              G4ForceCondition* condition);
       // It returns the MeanFreePath of the process for the current track :
       // (energy, material)
       // The previousStepSize and G4ForceCondition* are not used.
       // This function overloads a virtual function of the base class.		      
       // It is invoked by the ProcessManager of the Particle.
       
     G4double GetMicroscopicCrossSection(G4DynamicParticle* aDynamicGamma,
                                         G4Element*         anElement);
       // It returns the total CrossSectionPerAtom of the process, 
       // for the current DynamicGamma (energy), in anElement.					 
					  
     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step& aStep);
       // It computes the final state of the process (at end of step),
       // returned as a ParticleChange object.			    
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.
               
  protected:

     virtual G4double ComputeCrossSectionPerAtom(G4double GammaEnergy, 
                                                 G4double AtomicNumber);

     virtual G4double ComputeMeanFreePath(G4double GammaEnergy, 
                                          G4Material* aMaterial);
  private:
  
     // hide assignment operator as private 
     G4ComptonScattering& operator=(const G4ComptonScattering &right);
     G4ComptonScattering(const G4ComptonScattering& );
                                          
  private:
     
     G4PhysicsTable* theCrossSectionTable;    // table for crosssection
     G4PhysicsTable* theMeanFreePathTable;    // table for mean free path
       
     G4double LowestEnergyLimit ;      // low  energy limit of the crossection formula
     G4double HighestEnergyLimit ;     // high energy limit of the crossection formula
     G4int NumbBinTable ;              // number of bins in the crossection table
};

#include "G4ComptonScattering.icc"
  
#endif
 
