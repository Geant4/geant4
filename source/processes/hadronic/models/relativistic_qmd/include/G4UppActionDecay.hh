
#ifndef G4UPPACTIONDECAY_H
#define G4UPPACTIONDECAY_H


#include "G4VUppAction.hh"


class G4UppActionDecay : public G4VUppAction
{
public:

  G4UppActionDecay(const G4double time, 
		   G4UppTrack* pPtr);
  G4bool isValid() const;
  G4int Perform(const G4UppTrackVector& t) const { return 1; }
  G4int Perform(const G4UppTrackVector& t, G4UppInteraction& i) const;

private:

  G4UppTrack* particlePtr;

};


#endif // G4UPPACTIONDECAY_H
