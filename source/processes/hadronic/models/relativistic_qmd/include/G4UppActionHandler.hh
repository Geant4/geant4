//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#ifndef G4UPPACTIONHANDLER_H
#define G4UPPACTIONHANDLER_H


#include "G4UppInteraction.hh"
#include "G4VUppAction.hh"
#include "G4Scatterer.hh"
#include <g4std/vector>



class G4UppActionHandler
{
public:

  G4UppActionHandler(const G4UppTrackVector& allTracks);

  void addAction(G4VUppAction* anActionPtr);

  const G4VUppAction* getFirstAction() const;
  void deleteFirstAction();

  void updateActions();

  void cleanUp();

  G4bool empty() const 
    { return q.empty(); }

  void dump() const;

private:

  typedef G4VUppAction* valuetype;
  typedef vector<valuetype> queuetype;

  struct CmpAction {
    bool operator() (const G4VUppAction* a,const G4VUppAction* b) const;
  };

  queuetype q;
  const G4UppTrackVector* allTracksPtr;
  G4Scatterer aScatterer;


  G4double collisionTime(const G4UppTrack& i, 
			 const G4UppTrack& j) const;

  G4double minimumDistance(const G4UppTrack& i, 
			   const G4UppTrack& j) const;

  G4bool lastPartners(const G4UppTrack& i, 
		      const G4UppTrack& j) const;

  G4bool sameGroup(const G4UppTrack& i, 
		   const G4UppTrack& j) const;
};


#endif // G4UPPACTIONHANDLER_H
