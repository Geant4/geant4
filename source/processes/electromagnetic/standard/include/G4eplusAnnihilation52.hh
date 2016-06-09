//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4eplusAnnihilation52.hh,v 1.2 2006/06/29 19:52:18 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 10-01-97, crossection table + meanfreepath table, M.Maire
// 17-03-97, merge 'in fly' and 'at rest', M.Maire
// 31-08-98, new methods SetBining() and PrintInfo()
// 03-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 20-09-01, DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 05-08-04, suppress .icc file
// 13-08-04, public ComputeCrossSectionPerAtom() and ComputeMeanFreePath()   
// 09-11-04, Remove Store/Retrieve tables (V.Ivantchenko)
// 04-05-05, Add 52 to class name (V.Ivanchenko)
//

// class description
//
// e+ e- ---> gamma gamma
// inherit from G4VRestDiscreteProcess
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4eplusAnnihilation52_h
#define G4eplusAnnihilation52_h 1

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

class G4eplusAnnihilation52 : public G4VRestDiscreteProcess

{    
  public:  // with description
 
     G4eplusAnnihilation52(const G4String& processName ="annihil",
		               G4ProcessType type = fElectromagnetic);
 
    ~G4eplusAnnihilation52();

     G4bool IsApplicable(const G4ParticleDefinition&);
       // true for positron only.
            
     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);
       // Allows to define the binning of the PhysicsTables, 
       // before to build them.
            
     void BuildPhysicsTable(const G4ParticleDefinition&);
       // It builds the total CrossSectionPerAtom table, for e+,
       // and for every element contained in the elementTable.
       // It builds the MeanFreePath table, for e+,
       // and for every material contained in the materialTable.

     G4bool StorePhysicsTable(const G4ParticleDefinition* ,
			      const G4String& directory, G4bool);
       // store CrossSection and MeanFreePath tables into an external file
       // specified by 'directory' (must exist before invokation)

       //G4bool RetrievePhysicsTable(const G4ParticleDefinition* ,
       //				 const G4String& directory, G4bool);
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

     virtual
     G4double ComputeCrossSectionPerAtom(G4double PositKinEnergy,
                                         G4double AtomicNumber);

     G4double ComputeMeanFreePath(G4double PositKinEnergy, 
                                  G4Material* aMaterial);

  private:
  
   // hide assignment operator as private 
   G4eplusAnnihilation52& operator=(const G4eplusAnnihilation52& right);
   G4eplusAnnihilation52(const G4eplusAnnihilation52& );
      
  private:

     G4PhysicsTable* theCrossSectionTable;    
     G4PhysicsTable* theMeanFreePathTable;
     
     G4double LowestEnergyLimit;      // low  energy limit of the tables
     G4double HighestEnergyLimit;     // high energy limit of the tables 
     G4int NumbBinTable;              // number of bins in the tables
     G4double fminimalEnergy;         // minimalEnergy of produced particles     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif
 
