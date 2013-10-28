#ifndef Tst23ActionInit_h
#define Tst23ActionInit_h 1

#include "G4VUserActionInitialization.hh"

// forward declarations
//
class G4VUserPrimaryGeneratorAction;
class G4UserSteppingAction;

class Tst23ActionInit : public G4VUserActionInitialization
{

   public:
      
      // ctor & dtor
      Tst23ActionInit(); // : fPrimaryGen(0), fSteppingAct(0) {}
      virtual ~Tst23ActionInit() {}
      
      virtual void Build() const;
      virtual void BuildForMaster() const { return; }

      void SetAct( G4VUserPrimaryGeneratorAction* pg ){ fPrimaryGen=pg; return; }
      void SetAct( G4UserSteppingAction* sa )         { fSteppingAct=sa; return; }

   private:
            
      G4VUserPrimaryGeneratorAction* fPrimaryGen;
      G4UserSteppingAction*          fSteppingAct;
      
};

#endif
