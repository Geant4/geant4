
#include "Tst23ActionInit.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"

#include <assert.h>

Tst23ActionInit::Tst23ActionInit()
   : G4VUserActionInitialization(),
     fPrimaryGen(0),
     fSteppingAct(0)
{
}


void Tst23ActionInit::Build() const
{

   assert(fPrimaryGen);
   assert(fSteppingAct);
   
   SetUserAction( fPrimaryGen );
   SetUserAction( fSteppingAct );
   
   return;

}
