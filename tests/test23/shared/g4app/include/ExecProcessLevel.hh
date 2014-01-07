#ifndef ExecProcessLevel_H 
#define ExecProcessLevel_H 1

#include "ExecBase.hh"
#include "Beam.hh"
#include "TstTarget.hh"

// fwd declarations
// class TstTarget;
class G4Region;
class G4HadronicProcess;
class G4ProcessManager;
class G4Track;
class G4Step;
class G4StepPoint;

class ExecProcessLevel : public ExecBase
{

   public:
         
      ExecProcessLevel( const TstReader* pset ) : ExecBase(pset), fXSecOnTarget(0.) {}
      virtual ~ExecProcessLevel();
      
      virtual G4VParticleChange* DoEvent();
      
      const TstTarget* GetTarget() const { return fTarget; } 
      const Beam*      GetBeam()   const { return fBeam; }
      G4double         GetXSecOnTarget() const { return fXSecOnTarget; }
      
   protected:
   
      ExecProcessLevel();
      virtual void InitSetup(   const TstReader* );
      virtual void InitBeam(    const TstReader* );
      virtual void InitProcess( const TstReader* ) = 0;

      // protected data member
      G4HadronicProcess* fProcWrapper;
      G4double           fXSecOnTarget;
   
   private:
         
      // data members
      TstTarget*         fTarget;
      G4Region*          fRegion;
      G4ProcessManager*  fProcManager;
      Beam*              fBeam;
      G4Track*           fTrack;
      G4Step*            fStep;
      G4StepPoint*       fStepPoint;
      G4VParticleChange* fPartChange;

};

#endif
