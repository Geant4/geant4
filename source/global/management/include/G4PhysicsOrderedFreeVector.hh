// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsOrderedFreeVector.hh,v 1.1 1999-01-07 16:09:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// PhysicsOrderedFreeVector Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:	G4PhysicsOrderedFreeVector.hh
// Version:	1.0
// Created:     1996-08-13
// Author:      Juliet Armstrong
// Updated:     1997-03-25 by Peter Gumplinger
//		> cosmetics (only)
// mail:        gum@triumf.ca
//
// Description:
//    A physics ordered free vector inherits from G4PhysicsVector which
//    has values of energy-loss, cross-section, and other physics values
//    of a particle in matter in a given range of the energy, momentum,
//    etc.). In addition, the ordered free vector provides a method for
//    the user to insert energy/value pairs in sequence.  Methods to
//    Retrieve the Max and Min energies and values from the vector are
//    also provided. 
//
////////////////////////////////////////////////////////////////////////

#ifndef G4PhysicsOrderedFreeVector_h
#define G4PhysicsOrderedFreeVector_h 1

/////////////
// Includes
/////////////

#include <rw/tpordvec.h>
#include "G4PhysicsVector.hh"

/////////////////////
// Class Definition
/////////////////////

class G4PhysicsOrderedFreeVector : public G4PhysicsVector 
{
public:
	
        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

	G4PhysicsOrderedFreeVector();
	G4PhysicsOrderedFreeVector(G4double* Energies,
				   G4double* Values,
				   size_t VectorLength);

	~G4PhysicsOrderedFreeVector();

        ////////////
        // Methods
        ////////////

	void InsertValues(G4double energy, G4double value); 

	G4double GetLowEdgeEnergy(size_t binNumber) const;

	G4double GetMaxValue();

	G4double GetMinValue();

	G4double GetEnergy(G4double aValue);

	G4double GetMaxLowEdgeEnergy();

	G4double GetMinLowEdgeEnergy();

	void DumpValues();

private:

	size_t FindBinLocation(G4double theEnergy) const;

	size_t FindValueBinLocation(G4double aValue);

        G4double LinearInterpolationOfEnergy(G4double aValue, size_t theLocBin);
};

////////////////////
// Inline methods
////////////////////

inline
G4double G4PhysicsOrderedFreeVector::GetMaxValue()
{
	return dataVector.last();
}

inline
G4double G4PhysicsOrderedFreeVector::GetMinValue()
{
	return dataVector.first();
}

inline
G4double G4PhysicsOrderedFreeVector::GetMaxLowEdgeEnergy()
{
	return binVector.last();
}

inline
G4double G4PhysicsOrderedFreeVector::GetMinLowEdgeEnergy()
{
	return binVector.first();
}

inline
void G4PhysicsOrderedFreeVector::DumpValues()
{
   for (G4int i = 0; i < numberOfBin; i++) {
      G4cout << binVector[i] << "\t" << dataVector[i] << endl;
   }

}

inline
size_t G4PhysicsOrderedFreeVector::FindBinLocation(G4double theEnergy) const
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
   return (size_t)n1;
}

#endif /* G4PhysicsOrderedFreeVector_h */
