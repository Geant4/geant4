// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LPhysicsFreeVector.hh,v 1.1 1999-01-07 16:09:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------------
//
// Class G4LPhysicsFreeVector -- header file
//
// Derived from base class G4PhysicsVector
// This is a free vector for Low Energy Physics cross section data.
// The class name includes an "L" to distinguish it from other groups
// who may wish to implement a free vector in a different way.
// A subdivision method is used to find the energy|momentum bin.
//
// F.W. Jones, TRIUMF, 04-JUN-96
//
// 10-JUL-96 FWJ: adapted to changes in G4PhysicsVector.
//
// 27-MAR-97 FWJ: first version for Alpha release
// 20-JUN-97 FWJ: added comment re GetValue(): no longer virtual
//

#ifndef G4LPhysicsFreeVector_h
#define G4LPhysicsFreeVector_h 1

#include "G4PhysicsVector.hh"

class G4LPhysicsFreeVector : public G4PhysicsVector  
{
public:

   G4LPhysicsFreeVector();
   G4LPhysicsFreeVector(size_t nbin, G4double binmin, G4double binmax);

   ~G4LPhysicsFreeVector();

// G4PhysicsVector has PutValue() but it is inconvenient.
// Want to simultaneously fill the bin and data vectors.
   inline
   void PutValues(size_t binNumber, G4double binValue, G4double dataValue)
   {
      binVector(binNumber) = binValue;
      dataVector(binNumber) = dataValue;
   }

// Note that theEnergy could be energy, momentum, or whatever.
   inline
   G4double GetValue(G4double theEnergy, G4bool& isOutRange);

   inline
   void SetVerboseLevel(G4int value)
   {
      verboseLevel = value;
   }

   inline
   G4int GetVerboseLevel(G4int)
   {
      return verboseLevel;
   }

   inline
   G4double GetLastEnergy()
   {
      return lastEnergy;
   }

   inline
   size_t GetLastBin()
   {
      return lastBin;
   }

   void DumpValues();

private:

   G4int verboseLevel;

// Pure virtual in G4PhysicsVector
   inline
   size_t FindBinLocation(G4double theEnergy) const;
};


// Note: GetValue() is no longer virtual in the parent class
// G4PhysicsVector, so at present the following function cannot
// be called through a base class pointer.

inline 
G4double 
G4LPhysicsFreeVector::GetValue(G4double theEnergy, G4bool& isOutRange)
{
   G4double returnValue;

   //   verboseLevel = 2;

   if (theEnergy < edgeMin) {
      isOutRange = true;

      if (verboseLevel > 1) G4cout << "G4LPhysicsFreeVector::GetValue " <<
         theEnergy << " " << dataVector(0) << " " << isOutRange << endl;

      returnValue = dataVector(0);
   } 
   else if (theEnergy > edgeMax) {
      isOutRange = true;

      if (verboseLevel > 1) G4cout << "G4LPhysicsFreeVector::GetValue " <<
         theEnergy << " " << dataVector(numberOfBin - 1) << " " << 
         isOutRange << endl;

      returnValue = dataVector(numberOfBin - 1);
   } 
   else {
     isOutRange = false;

     G4int n = FindBinLocation(theEnergy);

     G4double dsde = (dataVector(n + 1) - dataVector(n))/
                     (binVector(n + 1) - binVector(n));

     if (verboseLevel > 1) G4cout << "G4LPhysicsFreeVector::GetValue " << 
        theEnergy << " " << dataVector(n) + (theEnergy - binVector(n))*dsde <<
        " " << isOutRange << endl;

     returnValue = dataVector(n) + (theEnergy - binVector(n))*dsde;
   }
   return returnValue;
}                                                                  


inline 
size_t
G4LPhysicsFreeVector::FindBinLocation(G4double theEnergy) const
{
   G4int n1 = 0;
   G4int n2 = numberOfBin/2;
   G4int n3 = numberOfBin - 1;
   while (n1 != n3 - 1) {
      if (theEnergy > binVector(n2))
         n1 = n2;
      else
         n3 = n2;
      n2 = n1 + (n3 - n1 + 1)/2;
   }
   if (verboseLevel > 1) G4cout << 
      "G4LPhysicsFreeVector::FindBinLocation:  returning " << n1 << endl;
   return (size_t)n1;
}
#endif
