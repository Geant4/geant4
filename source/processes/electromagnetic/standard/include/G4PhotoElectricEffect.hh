// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhotoElectricEffect.hh,v 1.5 1999-12-17 18:26:12 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ------------ G4PhotoElectricEffect physics process ------
//                   by Michel Maire, April 1996
//
// 12-06-96, Added SelectRandomAtom() method and new data member
//           for cumulative total cross section, by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 17-09-96, Dynamic array PartialSumSigma
//           split ComputeBindingEnergy(), M.Maire
// 08-01-97, crossection table + meanfreepath table, M.Maire
// 13-03-97, adapted for the new physics scheme, M.Maire
// 13-08-98, new methods SetBining() PrintInfo()
// 17-11-98, use table of atomic shells in PostStepDoIt, mma
// 06-01-99, Sandia crossSection below 50 keV, V.Grichine mma 
// ------------------------------------------------------------

// class description
//
// inherit from G4VDiscreteProcess
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4PhotoElectricEffect_h
#define G4PhotoElectricEffect_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ElementTable.hh"
#include "G4Gamma.hh" 
#include "G4Electron.hh"
#include "G4Step.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
class G4PhotoElectricEffect : public G4VDiscreteProcess
 
{
  public:  // with description
 
     G4PhotoElectricEffect(const G4String& processName ="phot");
 
    ~G4PhotoElectricEffect();

     G4bool IsApplicable(const G4ParticleDefinition&);
      // true for Gamma only.
           
     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);
       // Allows to define the binning of the PhysicsTables, 
       // before to build them.

     void BuildPhysicsTable(const G4ParticleDefinition& PhotonType);
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
     
     G4double GetMeanFreePath(const G4Track&          aTrack,
                                    G4double          previousStepSize,
                                    G4ForceCondition* condition);
       // It returns the MeanFreePath of the process for the current track :
       // (energy, material)
       // The previousStepSize and G4ForceCondition* are not used.
       // This function overloads a virtual function of the base class.		      
       // It is invoked by the ProcessManager of the Particle.
   
     G4double GetCrossSectionPerAtom(const G4DynamicParticle* aDynamicPhoton,
                                               G4Element*         anElement);
       // It returns the total CrossSectionPerAtom of the process, 
       // for the current DynamicGamma (energy), in anElement.
       
     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                     const G4Step&  aStep);
       // It computes the final state of the process (at end of step),
       // returned as a ParticleChange object.			    
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.
        
  protected:
  
     virtual inline G4double ComputeKBindingEnergy(G4double AtomicNumber);

     virtual inline G4double ComputeL1BindingEnergy(G4double AtomicNumber);

     virtual inline G4double ComputeL2BindingEnergy(G4double AtomicNumber);

     virtual G4double ComputeCrossSectionPerAtom(G4double PhotonEnergy, 
                                                 G4double AtomicNumber);
						 
     virtual G4double ComputeMeanFreePath(G4double PhotonEnergy, 
                                          G4Material* aMaterial);
					  
     G4double ComputeSandiaCrossSection (G4double PhotonEnergy, 
                                         G4double AtomicNumber);
						 					  
     inline G4double ComputeSandiaMeanFreePath(G4double PhotonEnergy, 
                                               G4Material* aMaterial);
					  
  private:

     G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicPhoton,
                                 G4Material* aMaterial);
  private:
  
     // hide assignment operator as private 
     G4PhotoElectricEffect& operator=(const G4PhotoElectricEffect &right);
     G4PhotoElectricEffect(const G4PhotoElectricEffect& );
       
  private:
  
     G4PhysicsTable* theCrossSectionTable;    // table for crossection
     G4PhysicsTable* theMeanFreePathTable;    // table for Mean free path
     
     G4double LowestEnergyLimit ;      // low  energy limit of the physics tables
     G4double HighestEnergyLimit ;     // high energy limit of the physics tables 
     G4int NumbBinTable ;              // number of bins in the tables

     G4double MeanFreePath;            // actual Mean Free Path (current medium)
};

#include "G4PhotoElectricEffect.icc"
  
#endif
 
