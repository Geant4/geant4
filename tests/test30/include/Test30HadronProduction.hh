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
// -------------------------------------------------------------
//      GEANT 4 class 
//
//      ---------- Test30HadronProduction -------
//                by Vladimir Ivanchenko, 12 March 2002 
// 
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#ifndef Test30HadronProduction_Test30HadronProduction_h
#define Test30HadronProduction_Test30HadronProduction_h 1


#include "Test30VSecondaryGenerator.hh"

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Step;
class G4Track;
class G4VParticleChange;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test30HadronProduction : public G4VDiscreteProcess
{
  public:     

     Test30HadronProduction(const G4String& processName = "HadronProduction" );

     ~Test30HadronProduction();

     void SetSecondaryGenerator(Test30VSecondaryGenerator*);

     G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                        			     G4double   previousStepSize,
			                             G4ForceCondition* condition);

     G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

     G4bool IsApplicable(const G4ParticleDefinition&) {return true;};
		 
		 G4double GetMass() {return theGenerator->GetMass();};


//  protected:

    
     G4double GetMeanFreePath(const G4Track& aTrack,
                             G4double   previousStepSize,
                             G4ForceCondition* condition) {return DBL_MAX;};
			    
  private:
  
     // hide assignment operator as private 
     Test30HadronProduction(const Test30HadronProduction&);
     Test30HadronProduction& operator = (const Test30HadronProduction &right);

     void InitializeMe();

     Test30VSecondaryGenerator* theGenerator;
     G4VParticleChange theChange;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif  // Test30HadronProduction_Test30HadronProduction_h

