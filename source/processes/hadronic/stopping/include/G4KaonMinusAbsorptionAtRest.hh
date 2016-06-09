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
//      File name:     G4KaonMinusAbsorptionAtRest.hh 
//
//      Author:        Christian V"olcker (Christian.Volcker@cern.ch),
// 
//      Creation date: 10. November 1997
//
// -------------------------------------------------------------------

#ifndef G4KaonMinusAbsorptionAtRest_h
#define G4KaonMinusAbsorptionAtRest_h 1

// Class Description:
//
// Process for nuclear absorption of K- at rest.
// To be used in your physics list in case you need this physics.

#include "globals.hh"
#include "Randomize.hh"
#include "G4VRestProcess.hh"
#include "G4ParticleTypes.hh"
#include "G4Nucleus.hh"
#include "G4DynamicParticle.hh"
#include "G4DynamicParticleVector.hh"
#include "G4NucleiProperties.hh"
#include "G4HadronicProcessType.hh"


// *********************************************************
class G4KaonMinusAbsorptionAtRest : public G4VRestProcess
// *********************************************************
{  
private:
  // hide assignment operator as private 
      G4KaonMinusAbsorptionAtRest& operator=(const G4KaonMinusAbsorptionAtRest &right);
      G4KaonMinusAbsorptionAtRest(const G4KaonMinusAbsorptionAtRest& );
public:
      G4KaonMinusAbsorptionAtRest(const G4String& processName ="KaonMinusAbsorptionAtRest", 
                       G4ProcessType   aType = fHadronic );
     ~G4KaonMinusAbsorptionAtRest();

  //override methods...
public: 
  G4bool IsApplicable(const G4ParticleDefinition& particle) {
               return( particle == *(G4KaonMinus::KaonMinus()) );
  }

  void PreparePhysicsTable(const G4ParticleDefinition&);

  void BuildPhysicsTable(const G4ParticleDefinition&);

  // the main method ...
     G4VParticleChange* AtRestDoIt(const G4Track& aTrack, const G4Step& aStep); 

protected:                         // why?? might be private....
  // zero mean lifetime
     G4double GetMeanLifeTime(const G4Track& aTrack,
			      G4ForceCondition* ) 
     {
     G4double result = 0;
     if(aTrack.GetMaterial()->GetNumberOfElements() == 1)
        if(aTrack.GetMaterial()->GetZ()<1.5) result = DBL_MAX;
     return result;
     }

private:
  // returns proton or neutron with fermi-momentum 
     G4DynamicParticle GetAbsorbingNucleon();

  // returns proton or neutron particle definition; 
     G4ParticleDefinition* SelectAbsorbingNucleon();
    
  // provides the neutron halo factor for absorption on nucleus surface. 
  // in the G4Nucleus
     G4double NeutronHaloFactor(G4double Z, G4double N);

  //  creates the reaction products
     G4DynamicParticleVector* KaonNucleonReaction();

  // secondary pion absorption in parent nucleus
  // if TRUE, then add excitation energy to the Nucleus
     G4bool AbsorbPionByNucleus(G4DynamicParticle* aPion);
     
  //  secondary Sigma-Lambda conversion
  // if conversion Done, then add excitation energy to the Nucleus
     G4DynamicParticle *SigmaLambdaConversion(G4DynamicParticle* aSigma);

  // instance variables ...
private:
  // pointer to current stopped hadron
     const G4DynamicParticle *stoppedHadron;

  // pointer to current target nucleus
     G4Nucleus* nucleus;

  // some constant parameters

     G4double pionAbsorptionRate;
     
  // primary production rates ( for absorption on Carbon)

     G4double rateLambdaZeroPiZero;
     G4double rateSigmaMinusPiPlus;
     G4double rateSigmaPlusPiMinus;
     G4double rateSigmaZeroPiZero;

     G4double rateLambdaZeroPiMinus;
     G4double rateSigmaZeroPiMinus;
     G4double rateSigmaMinusPiZero;


  // Sigma Lambda Conversion rates
  // for sigma- p -> lambda n
  //     sigma+ n -> lambda p
  //     sigma- n -> lambda 
  
     G4double sigmaPlusLambdaConversionRate; 
     G4double sigmaMinusLambdaConversionRate;
     G4double sigmaZeroLambdaConversionRate;
     
};

#endif

