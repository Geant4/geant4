
#ifndef G4VUPPACTION_H
#define G4VUPPACTION_H


#include "globals.hh"
#include "G4UppInteraction.hh"
#include "G4UppTrackVector.hh"


class G4VUppAction
{
public:

  virtual G4UppTrackChange* perform(const G4UppTrackVector& tracks) const = 0;

  virtual G4bool isValid() const = 0;

  virtual void setActionTime(const G4double newTime) 
    { actionTime = newTime; }
  virtual G4double getActionTime() const 
    { return actionTime; }

  virtual void dump() const 
    { G4cout << "Unknown Action" << G4endl; }

private:

  G4double actionTime;

};


#endif // G4VUPPACTION_H
