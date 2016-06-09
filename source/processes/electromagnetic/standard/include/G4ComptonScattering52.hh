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
// $Id: G4ComptonScattering52.hh,v 1.4 2007/05/16 14:00:56 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//------------------ G4ComptonScattering52 physics process -----------------------
//                   by Michel Maire, April 1996
//
// 10-06-96, updated by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 17-02-97, New Physics scheme
// 25-02-97, GetMeanFreePath() now is public function
// 12-03-97, new physics scheme again
// 13-08-98, new methods SetBining()  PrintInfo()
// 03-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 19-09-01, come back to previous ProcessName "compt"
// 20-09-01, DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 13-08-04, suppress icc file; make public ComputeCrossSectionPerAtom()   (mma)
// 09-11-04, Remove Retrieve tables (V.Ivantchenko)
// 04-05-05, Add 52 to class name (V.Ivanchenko)
// -----------------------------------------------------------------------------

// class description
//
// inherit from G4VDiscreteProcess
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4ComptonScattering52_h
#define G4ComptonScattering52_h 1

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

class G4ComptonScattering52 : public G4VDiscreteProcess

{
  public:  // with description

  G4ComptonScattering52(const G4String& processName ="compt",
		            G4ProcessType type = fElectromagnetic);

    ~G4ComptonScattering52();

     G4bool IsApplicable(const G4ParticleDefinition&);
       // true for Gamma only.

     void SetPhysicsTableBining(G4double lowE, G4double highE, G4int nBins);
       // Allows to define the binning of the PhysicsTables,
       // before to build them.

     void BuildPhysicsTable(const G4ParticleDefinition&);
       // It builds the total CrossSectionPerAtom table, for Gamma,
       // and for every element contained in the elementTable.
       // It builds the MeanFreePath table, for Gamma,
       // and for every material contained in the materialTable.

     G4bool StorePhysicsTable(const G4ParticleDefinition* ,
			      const G4String& directory, G4bool);
       // store CrossSection and MeanFreePath tables into an external file
       // specified by 'directory' (must exist before invokation)

       //G4bool RetrievePhysicsTable(const G4ParticleDefinition* ,
       //			 const G4String& directory, G4bool);
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

     G4double GetCrossSectionPerAtom(G4DynamicParticle* aDynamicGamma,
                                         G4Element*         anElement);
       // It returns the total CrossSectionPerAtom of the process, 
       // for the current DynamicGamma (energy), in anElement.					 
					  
     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step& aStep);
       // It computes the final state of the process (at end of step),
       // returned as a ParticleChange object.			    
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.

     virtual
     G4double ComputeCrossSectionPerAtom(G4double GammaEnergy, 
                                         G4double AtomicNumber);

     G4double ComputeMeanFreePath(G4double GammaEnergy, 
                                  G4Material* aMaterial);
  private:
  
     // hide assignment operator as private 
     G4ComptonScattering52& operator=(const G4ComptonScattering52 &right);
     G4ComptonScattering52(const G4ComptonScattering52& );
                                          
  private:
     
     G4PhysicsTable* theCrossSectionTable;    // table for crosssection
     G4PhysicsTable* theMeanFreePathTable;    // table for mean free path
       
     G4double LowestEnergyLimit;      // low  energy limit of the tables
     G4double HighestEnergyLimit;     // high energy limit of the tables
     G4int NumbBinTable;              // number of bins in the tables
     
  protected:   
     
     G4double fminimalEnergy;         // minimalEnergy of produced particles
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif
 
