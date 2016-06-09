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
// $Id: G4SynchrotronRadiation.hh,v 1.12 2004/11/10 08:53:18 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      
//      History: 
//      21-5-98  1 version , V. Grichine
//      28-05-01, V.Ivanchenko minor changes to provide ANSI -wall compilation 
// 
//
// ------------------------------------------------------------

#ifndef G4SynchrotronRadiation_h
#define G4SynchrotronRadiation_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4ThreeVector.hh"

#include "G4Track.hh"
#include "G4Step.hh"


#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"


#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"

class G4SynchrotronRadiation : public G4VDiscreteProcess
{
  public:

     G4SynchrotronRadiation(const G4String& processName =
			                    "SynchrotronRadiation",
		                  G4ProcessType type = fElectromagnetic);

    virtual ~G4SynchrotronRadiation();

  private:

     G4SynchrotronRadiation & operator=(const G4SynchrotronRadiation &right);

     G4SynchrotronRadiation(const G4SynchrotronRadiation&);

  public:  /////////////////    Post Step functions  //////////////////////////

     G4double GetMeanFreePath( const G4Track& track,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition ) ;

     G4VParticleChange *PostStepDoIt( const G4Track& track,
                                      const G4Step& Step    ) ;

     G4double GetPhotonEnergy( const G4Track& trackData,
                               const G4Step&  stepData      ) ;


     G4bool IsApplicable(const G4ParticleDefinition&);

     static G4double GetLambdaConst();
     static G4double GetEnergyConst();

  protected:

  private:

     static const G4double fLambdaConst ;

     static const G4double fEnergyConst ;

     static const G4double fIntegralProbabilityOfSR[200] ;


     const G4double
     LowestKineticEnergy;   // low  energy limit of the cross-section formula

     const G4double
     HighestKineticEnergy;  // high energy limit of the cross-section formula

     G4int TotBin;          // number of bins in the tables

     G4double CutInRange;

     const G4ParticleDefinition* theGamma; 
     const G4ParticleDefinition* theElectron;
     const G4ParticleDefinition* thePositron;

     const G4double* GammaCutInKineticEnergy;
     const G4double* ElectronCutInKineticEnergy;
     const G4double* PositronCutInKineticEnergy;
     const G4double* ParticleCutInKineticEnergy;


     G4double GammaCutInKineticEnergyNow;
     G4double ElectronCutInKineticEnergyNow;
     G4double PositronCutInKineticEnergyNow;
     G4double ParticleCutInKineticEnergyNow;
};

//////////////////////////  INLINE METHODS  /////////////////////////////

inline G4bool 
G4SynchrotronRadiation::IsApplicable( const G4ParticleDefinition& particle )
{
   return(   (&particle == (const G4ParticleDefinition *)theElectron)
           ||(&particle == (const G4ParticleDefinition *)thePositron)
         ) ;
}
  
#endif  // end of G4SynchrotronRadiation.hh
 
