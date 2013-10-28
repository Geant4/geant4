#ifndef ExecProcessLevel_H 
#define ExecProcessLevel_H 1

#include "ExecBase.hh"
#include "Beam.hh"

// fwd declarations
class TstTarget;
class G4Region;
class ProcessWrapper;
class G4ProcessManager;
class G4Track;
class G4Step;
class G4StepPoint;

class ExecProcessLevel : public ExecBase
{

   public:
      
      
      ExecProcessLevel( const TstReader* pset ) : ExecBase(pset) { InitSetup(pset); InitBeam(pset); }
      virtual ~ExecProcessLevel();
      
      virtual G4VParticleChange* DoEvent();
      
      const Beam* GetBeam() const { return fBeam; }
      
   protected:
   
      ExecProcessLevel();
      virtual void InitSetup( const TstReader* );
      virtual void InitBeam(  const TstReader* );

   private:
   
      void InitProcess( const TstReader* );
      
      // data members
      TstTarget*         fTarget;
      G4Region*          fRegion;
      ProcessWrapper*    fProcWrapper;
      G4ProcessManager*  fProcManager;
      Beam*              fBeam;
      G4Track*           fTrack;
      G4Step*            fStep;
      G4StepPoint*       fStepPoint;
      G4VParticleChange* fPartChange;

};

#endif
