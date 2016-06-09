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
// $Id: G4PhotoElectricEffect52.hh,v 1.1 2005/05/04 16:16:12 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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
// 03-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 19-09-01, come back to previous process name "phot"
// 20-09-01, DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 10-01-02, moved few function from icc to cc
// 17-04-02, Keep only Sandia crossSections. Remove BuildPhysicsTables.
//           Simplify public interface (mma)
// 29-04-02, Generate theta angle of the photoelectron from Sauter-Gavrila
//           distribution (mma)
// 13-08-04, suppress icc file; make public ComputeCrossSectionPerAtom()  (mma)
// 04-05-05, Add 52 to class name (V.Ivanchenko)
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// class description
//
// inherit from G4VDiscreteProcess
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4PhotoElectricEffect52_h
#define G4PhotoElectricEffect52_h 1

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
 
class G4PhotoElectricEffect52 : public G4VDiscreteProcess
 
{
  public:  // with description
 
     G4PhotoElectricEffect52(const G4String& processName ="phot",
		                 G4ProcessType type = fElectromagnetic);
 
    ~G4PhotoElectricEffect52();

     G4bool IsApplicable(const G4ParticleDefinition&);
      // true for Gamma only.
           
     void PrintInfoDefinition();
      // Print few lines of informations about the process.
     
     G4double GetMeanFreePath(const G4Track& aTrack,
                                    G4double previousStepSize,
                                    G4ForceCondition* condition);
       // It returns the MeanFreePath of the process for the current track :
       // (energy, material)
       // The previousStepSize and G4ForceCondition* are not used.
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.
       
   
     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                     const G4Step&  aStep);
       // It computes the final state of the process (at end of step),
       // returned as a ParticleChange object.			    
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.
        
  public:  // with description
  
       // utilities to access the cross Section and mean free path:
	   	  
     virtual
     G4double ComputeCrossSectionPerAtom(G4double PhotonEnergy, 
                                         G4double AtomicNumber);
		 
     virtual
     G4double ComputeMeanFreePath(G4double PhotonEnergy, 
                                  G4Material* aMaterial);
				  
  protected:
     virtual
     G4double ElecThetaDistribution(G4double ElecKineEnergy);				  

  private:

     G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicPhoton,
                                 G4Material* aMaterial);
  private:
  
     // hide assignment operator as private 
     G4PhotoElectricEffect52& operator=(const G4PhotoElectricEffect52 &right);
     G4PhotoElectricEffect52(const G4PhotoElectricEffect52& );
       
  private:
  
     G4double fminimalEnergy;      // minimalEnergy of produced particles
     
     G4double MeanFreePath;        // actual Mean Free Path (current medium)
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif
 
