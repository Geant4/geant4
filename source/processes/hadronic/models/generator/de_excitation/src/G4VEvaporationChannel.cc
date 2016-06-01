// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//

#include "G4VEvaporationChannel.hh"

G4VEvaporationChannel::G4VEvaporationChannel(const G4VEvaporationChannel &right)
{
 G4Exception("G4VEvaporationChannel::copy_constructor meant to not be accessable");
}




const G4VEvaporationChannel & G4VEvaporationChannel::operator=(const G4VEvaporationChannel &right)
{
  G4Exception("G4VEvaporationChannel::operator= meant to not be accessable");
  return *this;
}


G4bool G4VEvaporationChannel::operator==(const G4VEvaporationChannel &right) const
{
    return (this == (G4VEvaporationChannel *) &right);
    //  return false;
}

G4bool G4VEvaporationChannel::operator!=(const G4VEvaporationChannel &right) const
{
    return (this != (G4VEvaporationChannel *) &right);
    //  return true;
}


