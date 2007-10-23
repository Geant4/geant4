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
//      ---------- Test43HadronProduction -------
//           created from test30 files originally by Vladimir Ivanchenko
// 
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#ifndef Test43HadronProduction_Test43HadronProduction_h
#define Test43HadronProduction_Test43HadronProduction_h 1


#include "Test43VSecondaryGenerator.hh"

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4HadFinalState.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Step;
class G4Track;
class G4VParticleChange;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test43HadronProduction : public G4VDiscreteProcess
{
  public:

     Test43HadronProduction(const G4String& processName = "HadronProduction" );

     ~Test43HadronProduction();

     void SetSecondaryGenerator(Test43VSecondaryGenerator*);

     G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                        			         G4double   previousStepSize,
			                                 G4ForceCondition* condition);

     G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

     G4bool IsApplicable(const G4ParticleDefinition&) {return true;};

     G4double GetMass() {return theGenerator->GetMass();};


//  protected:


     G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*)
                                                        {return DBL_MAX;};

  private:

     // hide assignment operator as private
     Test43HadronProduction(const Test43HadronProduction&);
     Test43HadronProduction& operator = (const Test43HadronProduction &right);

     void InitializeMe();

     Test43VSecondaryGenerator* theGenerator;
     G4VParticleChange          theChange;

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif  
