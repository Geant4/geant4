//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
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
//              2000-11-11 by H.Kurashige
//              > use STL vector for dataVector and binVector
//              19 Jun. 2009-06-19 by V.Ivanchenko 
//              > removed hidden bin 
//
// mail:        gum@triumf.ca
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
  : G4PhysicsVector()
{
  type = T_G4PhysicsOrderedFreeVector;

  for (size_t i = 0 ; i < VectorLength ; i++)
    {
      InsertValues(Energies[i], Values[i]);
    }
}

G4PhysicsOrderedFreeVector::G4PhysicsOrderedFreeVector()
  : G4PhysicsVector()
{
  type = T_G4PhysicsOrderedFreeVector;
}

        ////////////////
        // Destructors
        ////////////////

G4PhysicsOrderedFreeVector::~G4PhysicsOrderedFreeVector() {}

        ////////////
        // Methods
        ////////////
  
void G4PhysicsOrderedFreeVector::InsertValues(G4double energy, G4double value)
{
        std::vector<G4double>::iterator binLoc =
                 std::lower_bound(binVector.begin(), binVector.end(), energy);

        size_t binIdx = binLoc - binVector.begin();	// Iterator difference!

        std::vector<G4double>::iterator dataLoc = dataVector.begin() + binIdx;

        binVector.insert(binLoc, energy);
        dataVector.insert(dataLoc, value);

        numberOfNodes++;
        edgeMin = binVector.front();
        edgeMax = binVector.back();
}

G4double  G4PhysicsOrderedFreeVector::GetLowEdgeEnergy(size_t binNumber) const
{
        return binVector[binNumber];
} 

G4double  G4PhysicsOrderedFreeVector::GetEnergy(G4double aValue)
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

size_t G4PhysicsOrderedFreeVector::FindValueBinLocation(G4double aValue)
{
   G4int n1 = 0;
   G4int n2 = numberOfNodes/2;
   G4int n3 = numberOfNodes - 1;
   while (n1 != n3 - 1) {
      if (aValue > dataVector[n2])
         { n1 = n2; }
      else
         { n3 = n2; }
      n2 = n1 + (n3 - n1 + 1)/2;
   }
   return (size_t)n1;
}

G4double G4PhysicsOrderedFreeVector::LinearInterpolationOfEnergy(G4double aValue,
								 size_t theLocBin)
{
  G4double intplFactor = (aValue-dataVector[theLocBin])
     / (dataVector[theLocBin+1]-dataVector[theLocBin]); // Interpolation factor

  return binVector[theLocBin] +
         ( binVector[theLocBin+1]-binVector[theLocBin] ) * intplFactor;
}


size_t G4PhysicsOrderedFreeVector::FindBinLocation(G4double theEnergy) const
{
   G4int n1 = 0;
   G4int n2 = numberOfNodes/2;
   G4int n3 = numberOfNodes - 1;
   while (n1 != n3 - 1)
   {
      if (theEnergy > binVector[n2])
         { n1 = n2; }
      else
         { n3 = n2; }
      n2 = n1 + (n3 - n1 + 1)/2;
   }
   return (size_t)n1;
}
