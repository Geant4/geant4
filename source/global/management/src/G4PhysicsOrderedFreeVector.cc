// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsOrderedFreeVector.cc,v 1.2 1999-11-16 17:46:52 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// PhysicsOrderedFreeVector Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4PhysicsOrderedFreeVector.cc
// Version:     2.0
// Created:     1996-08-13
// Author:      Juliet Armstrong
// Updated:     1997-03-25 by Peter Gumplinger
//              > cosmetics (only)
//              1998-11-11 by Peter Gumplinger
//              > initialize all data members of the base class in 
//                derived class constructors
// mail:        gum@triumf.ca
//
//
////////////////////////////////////////////////////////////////////////

#include "G4PhysicsOrderedFreeVector.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        /////////////////
        // Constructors
        /////////////////

G4PhysicsOrderedFreeVector::G4PhysicsOrderedFreeVector(G4double *Energies,
						       G4double *Values,
						       size_t VectorLength)
{
	
	ptrNextTable = 0;

        lastBin = INT_MAX;

	lastEnergy = -DBL_MAX;
	lastValue = DBL_MAX;

	numberOfBin = VectorLength;

	for (G4int i = 0 ; i < VectorLength ; i++)
	{
		binVector.insert(Energies[i]);	
		dataVector.insert(Values[i]); 
	}
        edgeMin = binVector.first();
        edgeMax = binVector.last();
}

G4PhysicsOrderedFreeVector::G4PhysicsOrderedFreeVector()
{
	ptrNextTable = 0;

        lastBin = INT_MAX;
        lastEnergy = -DBL_MAX;
        lastValue = DBL_MAX;

	edgeMin = 0.0;
	edgeMax = 0.0;
	numberOfBin = 0;
}

        ////////////////
        // Destructors
        ////////////////

G4PhysicsOrderedFreeVector::~G4PhysicsOrderedFreeVector() {}

        ////////////
        // Methods
        ////////////
  
void 
G4PhysicsOrderedFreeVector::InsertValues(G4double energy, G4double value)
{
        binVector.insert(energy);
        dataVector.insert(value);
	numberOfBin++;
        edgeMin = binVector.first();
        edgeMax = binVector.last();

}

G4double 
G4PhysicsOrderedFreeVector::GetLowEdgeEnergy(size_t binNumber) const
{
	return binVector[binNumber];
} 

G4double 
G4PhysicsOrderedFreeVector::GetEnergy(G4double aValue)
{

	if (aValue <= GetMinValue()) {
		return GetMinLowEdgeEnergy();
	} else if (aValue >= GetMaxValue()) {
		return GetMaxLowEdgeEnergy();
	} else { 
	size_t closestBin = FindValueBinLocation(aValue);
	G4double theEnergy = LinearInterpolationOfEnergy(aValue, closestBin);	

	return theEnergy;
	}
}

size_t
G4PhysicsOrderedFreeVector::FindValueBinLocation(G4double aValue)
{
   G4int n1 = 0;
   G4int n2 = numberOfBin/2;
   G4int n3 = numberOfBin - 1;
   while (n1 != n3 - 1) {
      if (aValue > dataVector(n2))
         n1 = n2;
      else
         n3 = n2;
      n2 = n1 + (n3 - n1 + 1)/2;
   }
   return (size_t)n1;
}

G4double 
G4PhysicsOrderedFreeVector::LinearInterpolationOfEnergy(G4double aValue,
						        size_t theLocBin)
{
  G4double intplFactor = (aValue-dataVector(theLocBin))
     / (dataVector(theLocBin+1)-dataVector(theLocBin)); // Interpolation factor

  return binVector(theLocBin) +
         ( binVector(theLocBin+1)-binVector(theLocBin) ) * intplFactor;
}
