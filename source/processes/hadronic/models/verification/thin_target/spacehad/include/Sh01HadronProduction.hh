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
//
//  
//
//  Sh01HadronProduction 
//             
// 
//    Modified: 
//
//    06.03.03  V. Grichine (based on V. Ivanchenko Test30)
//


#ifndef Sh01HadronProduction_Sh01HadronProduction_h
#define Sh01HadronProduction_Sh01HadronProduction_h 1


#include "Sh01SecondaryGenerator.hh"

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"



class G4Step;
class G4Track;
class G4VParticleChange;

//////////////////////////////////////////////////////////////////

class Sh01HadronProduction : public G4VDiscreteProcess
{
  public:     

     Sh01HadronProduction(const G4String& processName = "HadronProduction" );

     ~Sh01HadronProduction();

     void SetSecondaryGenerator(Sh01SecondaryGenerator*);

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
     Sh01HadronProduction(const Sh01HadronProduction&);
     Sh01HadronProduction& operator = (const Sh01HadronProduction &right);

     void InitializeMe();

     Sh01SecondaryGenerator* theGenerator;
     G4VParticleChange theChange;
};


///////////////////////////////////////////////////////////////////

#endif  // Sh01HadronProduction_Sh01HadronProduction_h

