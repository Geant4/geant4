
#ifndef G4VUPPACTION_H
#define G4VUPPACTION_H


#include "globals.hh"
#include "G4UppInteraction.hh"
#include "G4UppTrackVector.hh"


class G4VUppAction
{
public:

  virtual G4double getActionTime() const { return ActionTime; }
  virtual void setActionTime(const G4double time) { ActionTime = time; }
  virtual G4int Perform(const G4UppTrackVector& t) const = 0;
  virtual G4int Perform(const G4UppTrackVector& t, G4UppInteraction& i) const = 0;
  virtual G4bool isValid() const = 0;
  virtual void dump() const { G4cout << "Unknown Action" << G4endl; }
  virtual void dump(const G4UppTrackVector& t) const { G4cout << "Unknown Action" << G4endl; }

private:

  G4double ActionTime;

};


#endif // G4VUPPACTION_H
