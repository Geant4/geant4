
#ifndef G4VUPPFIELDTRANSPORT_H
#define G4VUPPFIELDTRANSPORT_H


#include "G4UppTrackVector.hh"


class G4VUppFieldtransport
{
public:

  virtual void Propagate(G4UppTrackVector& a, const G4double dTime) const = 0;

};


#endif // G4VUPPFIELDTRANSPORT_H
