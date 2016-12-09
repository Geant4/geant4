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
// $Id: G4PhysicsOrderedFreeVector.cc 98864 2016-08-15 11:53:26Z gcosmo $
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
//              2009-06-19 by V.Ivanchenko 
//              > removed hidden bin 
//              2013-10-02 by V.Ivanchenko removed FindBinLocation   
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4PhysicsOrderedFreeVector.hh"

G4PhysicsOrderedFreeVector::G4PhysicsOrderedFreeVector()
  : G4PhysicsVector()
{
        type = T_G4PhysicsOrderedFreeVector;
}

G4PhysicsOrderedFreeVector::G4PhysicsOrderedFreeVector(G4double *Energies,
                                                       G4double *Values,
                                                       size_t VectorLength)
  : G4PhysicsVector()
{
        type = T_G4PhysicsOrderedFreeVector;

        dataVector.reserve(VectorLength);
        binVector.reserve(VectorLength);

        for (size_t i = 0 ; i < VectorLength ; ++i)
        { 
          InsertValues(Energies[i], Values[i]);
        }
}

G4PhysicsOrderedFreeVector::~G4PhysicsOrderedFreeVector()
{}
  
void G4PhysicsOrderedFreeVector::InsertValues(G4double energy, G4double value)
{
        std::vector<G4double>::iterator binLoc =
                 std::lower_bound(binVector.begin(), binVector.end(), energy);

        size_t binIdx = binLoc - binVector.begin();	// Iterator difference!

        std::vector<G4double>::iterator dataLoc = dataVector.begin() + binIdx;

        binVector.insert(binLoc, energy);
        dataVector.insert(dataLoc, value);

        ++numberOfNodes;
        edgeMin = binVector.front();
        edgeMax = binVector.back();
}

G4double G4PhysicsOrderedFreeVector::GetEnergy(G4double aValue)
{
        G4double e;
        if (aValue <= GetMinValue()) {
          e = edgeMin;
        } else if (aValue >= GetMaxValue()) {
          e = edgeMax;
        } else { 
          size_t closestBin = FindValueBinLocation(aValue);
          e = LinearInterpolationOfEnergy(aValue, closestBin);
	}
        return e;
}

size_t G4PhysicsOrderedFreeVector::FindValueBinLocation(G4double aValue)
{
        size_t bin = std::lower_bound(dataVector.begin(), dataVector.end(), aValue)
                   - dataVector.begin() - 1;
        bin = std::min(bin, numberOfNodes-2);
        return bin;
}

G4double G4PhysicsOrderedFreeVector::LinearInterpolationOfEnergy(G4double aValue,
								 size_t bin)
{
        G4double res = binVector[bin];
        G4double del = dataVector[bin+1] - dataVector[bin];
        if(del > 0.0) { 
          res += (aValue - dataVector[bin])*(binVector[bin+1] - res)/del;  
        }
        return res;
}
