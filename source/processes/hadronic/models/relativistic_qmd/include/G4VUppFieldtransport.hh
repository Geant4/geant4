
#ifndef G4VUPPFIELDTRANSPORT_H
#define G4VUPPFIELDTRANSPORT_H


#include "G4UppTrackVector.hh"


class G4VUppFieldtransport
{
public:

  virtual void propagate(G4UppTrackVector& allTracks, 
			 const G4double dTime) const = 0;

  virtual string getName() const 
    { return "Unknown Fieldtransport"; }

};


#endif // G4VUPPFIELDTRANSPORT_H
