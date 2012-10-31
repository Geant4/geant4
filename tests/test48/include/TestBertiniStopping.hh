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
// -------------------------------------------------------------
//
//    Class: TestBertiniStopping
//    Purpose: attempt to interface Bertini Cascade for pi-, K- stopping/capture
//    Author: Julia Yarba
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#ifndef TestBertiniStopping_h
#define TestBertiniStopping_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestProcess.hh"
#include "G4HadFinalState.hh"
#include "G4HadronicProcessType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Step;
class G4Track;
class G4VParticleChange;
class G4Nucleus;
class G4HadronicInteraction;
class G4CascadeInterface;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class TestBertiniStopping : public G4VRestProcess
{
  public:

     TestBertiniStopping(const G4String&     processName = "HadronBertiniStopping",
     			       G4ProcessType processType = fHadronic );
     

     ~TestBertiniStopping();

     G4double AtRestGetPhysicalInteractionLength(const G4Track&,
						       G4ForceCondition*) { return 0.; } // since mean lifetime is 0.
     // zero mean lifetime
     G4double GetMeanLifeTime(const G4Track& ,
			            G4ForceCondition* ) {return 0.0;}

     G4VParticleChange* AtRestDoIt(const G4Track& aTrack, const G4Step& aStep);
     
     G4bool IsApplicable(const G4ParticleDefinition&) {return true;};

     // return number of secondaries produced
     G4int GetNumberOfSecondaries() { return fNSec; }

     void InitTarget( G4Material* );
     void UsePreCompound();
     
  private:

     // hide assignment operator as private
     TestBertiniStopping( const TestBertiniStopping& );
     TestBertiniStopping& operator = ( const TestBertiniStopping& right );
     
     // G4HadronicInteraction*     fModel;
     G4CascadeInterface*        fModel;
     G4Nucleus*                 fTargetNucleus;
     G4VParticleChange          fPartChange;     
     G4int                      fNSec;

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif  // TestBertiniStopping_h

