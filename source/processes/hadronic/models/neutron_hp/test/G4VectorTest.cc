// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VectorTest.cc,v 1.2 1999-12-15 14:53:22 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "g4templates.hh"
#include "G4NeutronHPVector.hh"
#include "../src/G4NeutronHPVector.cc"

main()
{

  G4NeutronHPVector aVector;
  
  aVector.SetData(0, .000000E+00, .000000E+00);
  aVector.SetData(1, .500000E+07, .000000E+00);
  aVector.SetData(2, .550000E+07, .465857E-07);
  aVector.SetData(3, .550000E+07, .465857E-07);
  aVector.SetData(4, .600000E+07, .727875E-07);
  aVector.SetData(5, .700000E+07, .107729E-06);
  aVector.SetData(6, .900000E+07, .116469E-06);
  aVector.SetData(7, .140000E+08, .110639E-06);
  aVector.SetData(8, .149000E+08, .582307E-07);
  aVector.SetData(9, .149100E+08, .000000E+00);

  for(G4int i=0; i<100000000000; i++)
  {
    G4double it = aVector.Sample();
    G4cout << it<<G4endl;
  }
}
