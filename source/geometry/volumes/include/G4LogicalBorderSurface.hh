// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalBorderSurface.hh,v 1.2 1999-11-11 15:35:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// G4LogicalBorderSurface Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4LogicalBorderSurface.hh
// Description: A Logical Surface class for surfaces defined by the
//               boundary of two physical volumes.
// Version:     1.0
// Created:     1997-06-17
// Author:      John Apostolakis
// mail:        John.Apostolakis@cern.ch
// Modified:    1997-06-16  John Apostolakis
//
// Id tag:      
////////////////////////////////////////////////////////////////////////

#ifndef G4LogicalBorderSurface_h
#define G4LogicalBorderSurface_h 1

/////////////
// Includes
/////////////

#include  "G4LogicalSurface.hh"
#include "G4VPhysicalVolume.hh"

// G4RWTPtrOrderedVector
#include "g4rw/tpordvec.h"

class G4Event;

class G4VPhysicalVolume;

/////////////////////
// Class Definition
/////////////////////

class G4LogicalBorderSurface: public G4LogicalSurface
{
        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////
public:

	G4LogicalBorderSurface(const G4String& name,
			       G4VPhysicalVolume* vol1, 
			       G4VPhysicalVolume* vol2,
			       G4OpticalSurface* opticsSurface);

	~G4LogicalBorderSurface();

	////////////
	// Methods
        ////////////
public:
	// public methods

	static G4LogicalBorderSurface* GetSurface(const G4VPhysicalVolume* vol1,
					    const G4VPhysicalVolume* vol2);

	void       SetPhysicalVolumes(G4VPhysicalVolume* vol1,
				      G4VPhysicalVolume* vol2)
				       { Volume1 = vol1; Volume2 = vol2; }

	G4VPhysicalVolume* GetVolume1() const {return Volume1;}
	G4VPhysicalVolume* GetVolume2() const {return Volume2;}

        // These are potentially dangerous.
	void		   SetVolume1(G4VPhysicalVolume* vol1)
						{Volume1 = vol1;}

	void		   SetVolume2(G4VPhysicalVolume* vol2)
						{Volume2 = vol2;}

	//   Methods dealing with the table of surfaces.
	//
        static const G4RWTPtrOrderedVector<G4LogicalBorderSurface>* GetSurfaceTable()
					   { return &theBorderSurfaceTable; }
        static size_t GetNumberOfBorderSurfaces()
				   { return theBorderSurfaceTable.length(); }
	static void DumpInfo(); 

	size_t GetIndex() const { return theIndexInTable; }

        //////////////
        // Operators
        //////////////
public:
	G4int operator==(const G4LogicalBorderSurface &right) const;
	G4int operator!=(const G4LogicalBorderSurface &right) const;

private:
        G4LogicalBorderSurface(const G4LogicalBorderSurface &right);
	const G4LogicalBorderSurface & operator=(const G4LogicalBorderSurface &right);

	// ------------------
	// Basic data members ( To define a 'logical' surface)
	// ------------------

private:
	G4VPhysicalVolume* Volume1;	// Physical Volume pointer on side 1
	G4VPhysicalVolume* Volume2;	// Physical Volume pointer on side 2

//	The static Table of Surfaces
	static G4RWTPtrOrderedVector<G4LogicalBorderSurface> theBorderSurfaceTable;

	size_t theIndexInTable;		// Index of surface in the surface table

};

typedef G4RWTPtrOrderedVector<G4LogicalBorderSurface> G4LogicalBorderSurfaceTable;

////////////////////
// Inline methods
////////////////////

#endif /* G4LogicalBorderSurface_h */
