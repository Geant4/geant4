//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4eplusAnnihilation.hh,v 1.6 2001-08-09 17:24:23 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 10-01-97, crossection table + meanfreepath table, M.Maire
// 17-03-97, merge 'in fly' and 'at rest', M.Maire
// 31-08-98, new methods SetBining() and PrintInfo()
// 03-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)   
// 

// class description
//
// e+ e- ---> gamma gamma
// inherit from G4VRestDiscreteProcess
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4eplusAnnihilation_h
#define G4eplusAnnihilation_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VRestDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh" 
#include "G4ElementTable.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4eplusAnnihilation : public G4VRestDiscreteProcess
 
{    
  public:  // with description
 
     G4eplusAnnihilation(const G4String& processName ="annihil");
 
    ~G4eplusAnnihilation();

     G4bool IsApplicable(const G4ParticleDefinition&);
       // true for positron only.
            
     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);
       // Allows to define the binning of the PhysicsTables, 
       // before to build them.
            
     void BuildThePhysicsTable();
       // It builds the total CrossSectionPerAtom table, for e+,
       // and for every element contained in the elementTable.
       // It builds the MeanFreePath table, for e+,
       // and for every material contained in the materialTable.       
       // It is invoked by the constructor. 
       
     G4bool StorePhysicsTable(G4ParticleDefinition* ,
			      const G4String& directory, G4bool);
       // store CrossSection and MeanFreePath tables into an external file
       // specified by 'directory' (must exist before invokation)

     G4bool RetrievePhysicsTable(G4ParticleDefinition* ,
				 const G4String& directory, G4bool);
       // retrieve CrossSection and MeanFreePath tables from an external file
       // specified by 'directory' 
                   
     void PrintInfoDefinition();
       // Print few lines of informations about the process: validity range,
       // origine ..etc..
       // Invoked by BuildThePhysicsTable(). 
           
     G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double previousStepSize,
                              G4ForceCondition* condition);
       // It returns the MeanFreePath of the process for the current track :
       // (energy, material)
       // The previousStepSize and G4ForceCondition* are not used.
       // This function overloads a virtual function of the base class.		      
       // It is invoked by the ProcessManager of the Particle.
        
     G4double GetCrossSectionPerAtom(G4DynamicParticle* aDynamicPositron,
                                     G4Element*         anElement);     
       // It returns the total CrossSectionPerAtom of the process, 
       // for the current DynamicPositron (energy), in anElement.

     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step& aStep); 
       // It computes the final state of the process (at end of step),
       // returned as a ParticleChange object.			    
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.

     G4double GetMeanLifeTime(const G4Track& aTrack,
                              G4ForceCondition* condition);
       // It is invoked by the ProcessManager of the Positron if this
       // e+ has a kinetic energy null. Then it return 0 to force the
       // call of AtRestDoIt.
       // This function overloads a virtual function of the base class.
              
     G4VParticleChange* AtRestDoIt(const G4Track& aTrack,
                                   const G4Step& aStep); 
       // It computes the final state of the process:
       //          e+ (at rest) e- (at rest)  ---> gamma gamma,
       // returned as a ParticleChange object.			    
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.

  protected:

     virtual G4double ComputeCrossSectionPerAtom(G4double PositKinEnergy,
                                                 G4double AtomicNumber);

     virtual G4double ComputeMeanFreePath(G4double PositKinEnergy, 
                                          G4Material* aMaterial);

  private:
  
   // hide assignment operator as private 
   G4eplusAnnihilation& operator=(const G4eplusAnnihilation& right);
   G4eplusAnnihilation(const G4eplusAnnihilation& );
      
  private:

     G4PhysicsTable* theCrossSectionTable;    
     G4PhysicsTable* theMeanFreePathTable;
     
     G4double LowestEnergyLimit;      // low  energy limit of the tables
     G4double HighestEnergyLimit;     // high energy limit of the tables 
     G4int NumbBinTable;              // number of bins in the tables
};

#include "G4eplusAnnihilation.icc"
  
#endif
 
