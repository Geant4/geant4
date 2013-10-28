#ifndef Tst23ParticleChange_H
#define Tst23ParticleChange_H 1

#include "G4VParticleChange.hh"

class Tst23ParticleChange : public G4VParticleChange
{

   public:
   
      Tst23ParticleChange() : fIsFirstInter(false) {}
      Tst23ParticleChange( bool isFirst ) : fIsFirstInter( isFirst ) {}
      virtual ~Tst23ParticleChange() {}
      
      bool IsFisrtInteraction() const { return fIsFirstInter; }

   private:
   
      bool fIsFirstInter;

};

#endif
