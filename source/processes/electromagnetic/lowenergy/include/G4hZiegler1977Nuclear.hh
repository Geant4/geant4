// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// For information related to this code contact:
// Geant4 Collaboration
//
// File name:     G4hZiegler1977Nuclear
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 20 July 2000
//
// Modifications: 
// 20/07/2000  V.Ivanchenko First implementation
//
// Class Description: 
//
// Nuclear stopping power parametrised according to
// J.F.Ziegler, Helium Stopping Powers and
// Ranges in All Elemental Matter, Vol.4, Pergamon Press, 1977
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#ifndef G4hZiegler1977Nuclear_h
#define G4hZiegler1977Nuclear_h 1

#include "G4VhNuclearStoppingPower.hh"

class G4hZiegler1977Nuclear : public G4VhNuclearStoppingPower
{

public:

  G4hZiegler1977Nuclear();

  ~G4hZiegler1977Nuclear();

  G4double NuclearStoppingPower(G4double kineticEnergy,
                                G4double z1, G4double z2, 
                                G4double m1, G4double m2) const;
 
protected:

private:

};

#endif
