// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LPhysicsFreeVector.cc,v 1.1 1999-01-07 16:09:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// Class G4LPhysicsFreeVector
// Derived from base class G4PhysicsVector
// This is a free vector for Low Energy Physics cross section data
//
// F.W. Jones, TRIUMF, 04-JUN-96
//
// 27-MAR-97 FWJ: first version for Alpha release
//

#include "G4LPhysicsFreeVector.hh"

#include <stdio.h>

G4LPhysicsFreeVector::G4LPhysicsFreeVector()
{
   ptrNextTable = 0;
   edgeMin = 0.0;
   edgeMax = 0.0;
   numberOfBin = 0;
   verboseLevel = 0;
}

G4LPhysicsFreeVector::G4LPhysicsFreeVector(size_t nbin, G4double binmin,
                                           G4double binmax)
{
   edgeMin = binmin;
   edgeMax = binmax;
   numberOfBin = nbin;
   lastEnergy = 0.;
   lastValue = 0.;
   lastBin = 0;
   binVector.resize(nbin);
   dataVector.resize(nbin);
   ptrNextTable = 0;
   verboseLevel = 0;
}  

G4LPhysicsFreeVector::~G4LPhysicsFreeVector()
{
}

void G4LPhysicsFreeVector::DumpValues()
{
   for (G4int i = 0; i < numberOfBin; i++) {
     //      printf(" %12.4f   %7.1f\n", binVector(i), dataVector(i)*1.e-27);
      printf(" %12.4f   %7.1f\n", binVector(i), dataVector(i)/millibarn);
   }
}
