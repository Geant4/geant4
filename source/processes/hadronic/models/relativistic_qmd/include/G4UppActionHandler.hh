
#ifndef G4UPPACTIONHANDLER_H
#define G4UPPACTIONHANDLER_H


#include "G4UppInteraction.hh"
#include "G4VUppAction.hh"
#include <g4std/vector>


struct CmpAction {
  bool operator() (const G4VUppAction* a,const G4VUppAction* b) const;
};


class G4UppActionHandler
{
public:

  const G4VUppAction& getFirstAction() const;
  void deleteFirstAction();
  void UpdateActions(const G4UppTrackVector& v);
  void CleanUp();
  void addAction(G4VUppAction* a);
  G4bool empty() const { return q.empty(); }
  void dump(const G4UppTrackVector& t) const;

private:

  typedef G4VUppAction* valuetype;
  typedef vector<valuetype> queuetype;

  queuetype q;

  G4double CollisionTime(const G4UppTrack& i, 
			 const G4UppTrack& j) const;
  G4double MinimumDistance(const G4UppTrack& i, 
			   const G4UppTrack& j) const;
  G4bool lastPartners(const G4UppTrack& i, 
		      const G4UppTrack& j) const;
  G4bool sameGroup(const G4UppTrack& i, 
		   const G4UppTrack& j) const;
};


#endif // G4UPPACTIONHANDLER_H
