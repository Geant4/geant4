// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserSpecialCuts.hh,v 1.3 1999-12-15 14:53:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//                  15 April 1998  M.Maire
// ------------------------------------------------------------
//
// class description
//
// a "process" to be registered to the process manager of each particle,
// in the UserPhysicsList, in order to take account the user limits defined in the 
// G4UserLimits attached to a logical volume.
//

#ifndef G4UserSpecialCuts_h
#define G4UserSpecialCuts_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4UserSpecialCuts : public G4VProcess 
{
  public:  //with description     

     G4UserSpecialCuts(const G4String& processName ="UserSpecialCut" );

    ~G4UserSpecialCuts();

     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

     virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );
			    
  public:  //without description 
     			    
     //  no operation in  AtRestGPIL
     virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4ForceCondition* 
			    ){ return -1.0; };
			    
     //  no operation in  AtRestDoIt      
     virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    ){return NULL;};

     //  no operation in  AlongStepGPIL
     virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
			     G4double  ,
			     G4double  ,
			     G4double& ,
                             G4GPILSelection*
			    ){ return -1.0; };

     //  no operation in  AlongStepDoIt
     virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    ) {return NULL;};

  private:
  
  // hide assignment operator as private 
      G4UserSpecialCuts(G4UserSpecialCuts&);
      G4UserSpecialCuts& operator=(const G4UserSpecialCuts& right);

};

#endif

