// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LPhysicsFreeVector.hh,v 1.6 2001-03-09 12:08:18 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------------
//
// Class G4LPhysicsFreeVector -- header file
//
// Class description:
//
// Derived from base class G4PhysicsVector
// This is a free vector for Low Energy Physics cross section data.
// The class name includes an "L" to distinguish it from other groups
// who may wish to implement a free vector in a different way.
// A subdivision method is used to find the energy|momentum bin.

// F.W. Jones, TRIUMF, 04-JUN-96
//
// 10-JUL-96 FWJ: adapted to changes in G4PhysicsVector.
//
// 27-MAR-97 FWJ: first version for Alpha release
// 20-JUN-97 FWJ: added comment re GetValue(): no longer virtual
// 11-NOV-00 H.Kurashige: use STL vector for dataVector and binVector
//

#ifndef G4LPhysicsFreeVector_h
#define G4LPhysicsFreeVector_h 1

#include "G4PhysicsVector.hh"

class G4LPhysicsFreeVector : public G4PhysicsVector  
{

public: 

   G4LPhysicsFreeVector();

public: // with description

   G4LPhysicsFreeVector(size_t nbin, G4double binmin, G4double binmax);

   ~G4LPhysicsFreeVector();

   void PutValues(size_t binNumber, G4double binValue, G4double dataValue);
     // G4PhysicsVector has PutValue() but it is inconvenient.
     // Want to simultaneously fill the bin and data vectors.

   G4double GetValue(G4double theEnergy, G4bool& isOutRange);
     // Note that theEnergy could be energy, momentum, or whatever.

   void SetVerboseLevel(G4int value);

   G4int GetVerboseLevel(G4int);

   G4double GetLastEnergy();

   size_t GetLastBin();

   void DumpValues();

private:

   G4int verboseLevel;

   size_t FindBinLocation(G4double theEnergy) const;
     // Pure virtual in G4PhysicsVector
};

#include "G4LPhysicsFreeVector.icc"

#endif
