#ifndef Tst23ParticleChange_H
#define Tst23ParticleChange_H 1

#include "G4VParticleChange.hh"

class Tst23ParticleChange : public G4VParticleChange
{

   public:
   
      Tst23ParticleChange() : fIsFirstInter(false), fIncomingTrack(0) {}
      Tst23ParticleChange( bool isFirst ) : fIsFirstInter( isFirst ), fIncomingTrack(0) {}
      virtual ~Tst23ParticleChange() {}
      
      void           SetIncomingTrack( G4Track* trk ) { fIncomingTrack=trk; return; }
      bool           IsFisrtInteraction() const { return fIsFirstInter; }
      const G4Track* GetIncomingTrack()   const { return fIncomingTrack; }

   private:
   
      bool     fIsFirstInter;
      G4Track* fIncomingTrack;

};

#endif
