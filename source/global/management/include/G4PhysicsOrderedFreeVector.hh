//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PhysicsOrderedFreeVector.hh,v 1.8 2001-07-11 10:00:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// PhysicsOrderedFreeVector Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:	G4PhysicsOrderedFreeVector.hh
// Created:     1996-08-13
// Author:      Juliet Armstrong
// Updated:     1997-03-25 by Peter Gumplinger
//		> cosmetics (only)
//              2000-11-11 by H.Kurashige
//              > use STL vector for dataVector and binVector
// mail:        gum@triumf.ca
//
// Class description:
//
//    A physics ordered free vector inherits from G4PhysicsVector which
//    has values of energy-loss, cross-section, and other physics values
//    of a particle in matter in a given range of the energy, momentum,
//    etc.). In addition, the ordered free vector provides a method for
//    the user to insert energy/value pairs in sequence.  Methods to
//    Retrieve the Max and Min energies and values from the vector are
//    also provided. 

////////////////////////////////////////////////////////////////////////

#ifndef G4PhysicsOrderedFreeVector_h
#define G4PhysicsOrderedFreeVector_h 1

/////////////
// Includes
/////////////

#include "G4PhysicsVector.hh"

/////////////////////
// Class Definition
/////////////////////

class G4PhysicsOrderedFreeVector : public G4PhysicsVector 
{

  public: // with description
	
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

#include "G4PhysicsOrderedFreeVector.icc"

#endif /* G4PhysicsOrderedFreeVector_h */






