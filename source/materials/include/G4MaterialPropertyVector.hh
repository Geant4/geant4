// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MaterialPropertyVector.hh,v 1.2 1999-10-30 01:42:31 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// G4MaterialPropertyVector Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4MaterialPropertyVector.hh
//
// Description: A one-to-one mapping from Photon Momentum to some
//              optical property 
// Version:     1.0
// Created:     1996-02-08
// Author:      Juliet Armstrong
// Updated:     1999-10-29 add method and class decriptors
//              1997-03-25 by Peter Gumplinger
//              > value.h -> templates.hh
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4MaterialPropertyVector_h   
#define G4MaterialPropertyVector_h 1

/////////////
// Includes
/////////////

#include <rw/tpsrtvec.h>
#include <rw/cstring.h>
#include "G4MPVEntry.hh"

// Class Description:
// A one-to-one mapping from Photon Momentum to some optical property.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4MaterialPropertyVector {

public: // Without description

	//////////////
	// Operators
	//////////////

	G4bool operator ++();
        G4MaterialPropertyVector&
                   operator =(const G4MaterialPropertyVector &right);

	/////////////////
	// Constructors
	/////////////////

	G4MaterialPropertyVector() : MPV(0) 
	{
		CurrentEntry = -1;
		NumEntries   = 0;
	};

public: // With description
	
	G4MaterialPropertyVector(G4double *PhotonMomenta, 
		  	   	 G4double *PropertyValues,
				 G4int     NumElements);
        // Constructor of G4MaterialPropertyVector object.

public: // Without description

	G4MaterialPropertyVector(const G4MaterialPropertyVector &right);

        ///////////////
        // Destructor
        ///////////////

	~G4MaterialPropertyVector();

        ////////////
        // Methods
        ////////////

public: // With description

	void ResetIterator();

        void AddElement(G4double aPhotonMomentum, 
			G4double aPropertyValue);
        // Add a new element (pair of numbers) to the G4MaterialPropertyVector.
	void RemoveElement(G4double aPhotonMomentum);
        // Remove the element with given x-value.

        G4double GetProperty(G4double aPhotonMomentum) const;
        // Returns the y-value for given x-value (with interpolation).
	G4double GetPhotonMomentum(G4double aProperty) const;
        // Returns the x-value for given y-value (with interpolation).
        // NOTE: Assumes that the y-value is an increasing function of 
        //       the x-value. Returns the x-value corresponding to the 
        //       y-value passed in. If several x-values correspond to 
        //       the y-value passed in, the method returns the first 
        //       x-value in the vector that corresponds to that value.
	// For use with G4MaterialPropertyVector iterator: return
        // property (or Photon momentum) at current point of iterator.

	G4double GetProperty() const;
	G4double GetPhotonMomentum() const;

	G4double GetMaxProperty() const;
	G4double GetMinProperty() const;
	G4double GetMaxPhotonMomentum() const;
	G4double GetMinPhotonMomentum() const;
		
	//////////
	// Tests
	//////////

	void DumpVector();	

private:

	/////////////////////
	// Helper Functions
	/////////////////////

	G4MPVEntry GetEntry(G4int i) const;

	void GetAdjacentBins(G4double aPhotonMomentum,
                             G4int *left,G4int *right) const;

	/////////////////////////
        // Private Data Members
	/////////////////////////

	RWTPtrSortedVector<G4MPVEntry> MPV;
	G4int NumEntries;
	G4int CurrentEntry;
};

///////////////////
// Inline methods
///////////////////

inline 
G4double G4MaterialPropertyVector::GetMaxProperty() const
{
	return MPV.last()->GetProperty();
}

inline
G4double G4MaterialPropertyVector::GetMinProperty() const
{
	return MPV.first()->GetProperty();
}

inline
G4double G4MaterialPropertyVector::GetMaxPhotonMomentum() const
{
	return MPV.last()->GetPhotonMomentum();
}

inline 
G4double G4MaterialPropertyVector::GetMinPhotonMomentum() const
{
	return MPV.first()->GetPhotonMomentum();
}

#endif /* G4MaterialPropertyVector_h */
