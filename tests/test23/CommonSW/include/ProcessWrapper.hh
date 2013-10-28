#ifndef ProcessWrapper_HH
#define ProcessWrapper_HH 1

#include "G4ios.hh"
#include "globals.hh"

#include "G4VDiscreteProcess.hh"

#include "G4HadronicInteraction.hh"
#include "G4VPartonStringModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4ExcitedStringDecay.hh"

// forward declarations
//
//class G4Step;
//class G4Track;

// NOTE: In principle, there's a class G4WrapperProcess ( : public G4VProcess)
// under /processes/management, but that one seems to be aiming something different'
// this one is to interface High Energy models/processes

class ProcessWrapper : public G4VDiscreteProcess
{

   public:
   
     // ctor & dtor
     ProcessWrapper( const G4String& name = "ProcessWrapper",  
                           G4ProcessType processType = fHadronic ) : G4VDiscreteProcess(name,processType), 
                                                                     fInteractionModel(0), fStringModel(0),
                                                                     fCascade(0), fStringDecay(0), fUseLundStrFragm(false) {}
     virtual ~ProcessWrapper();
     
     void UseG4LundStringFragm( bool g4lund=true ) { fUseLundStrFragm=g4lund; return ; } 
     virtual void Compose() = 0;
     
     G4double PostStepGetPhysicalInteractionLength(const G4Track&,
                        			         G4double,
			                                 G4ForceCondition* condition);

     virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

     G4bool IsApplicable(const G4ParticleDefinition&) {return true;};

     // G4double GetMass() {return theGenerator->GetMass();};


     G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*)
                                                        {return DBL_MAX;};

     void InitTarget( G4Material* mat ) { fTargetNucleus = new G4Nucleus(mat); return; }
  
   protected:
      
      // in principle, I can make it directly G4TheoFSGenerator*
      // because it's the same among FTF(p), QGS(P), and QGS(B)...
      //
      G4HadronicInteraction*           fInteractionModel;
      G4VPartonStringModel*            fStringModel;
      G4GeneratorPrecompoundInterface* fCascade; // need to check what type is G4BinaryCascade &
      G4ExcitedStringDecay*            fStringDecay;
      bool                             fUseLundStrFragm;
      
      G4Nucleus*                       fTargetNucleus;
      
      G4VParticleChange                fPartChange;
                                                 
};

#endif
