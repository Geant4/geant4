// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalBorderSurface.hh,v 1.5 2000-11-01 16:51:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// class G4LogicalBorderSurface
////////////////////////////////////////////////////////////////////////
//
// Class description:
//
// A Logical Surface class for surfaces defined by the boundary
// of two physical volumes.

// History:
// -------
// Created:     1997-06-17
// Author:      John Apostolakis
// mail:        John.Apostolakis@cern.ch
// Modified:    1997-06-16  John Apostolakis
//
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

	static G4LogicalBorderSurface* GetSurface(const G4VPhysicalVolume* vol1,
					          const G4VPhysicalVolume* vol2);

	inline void  SetPhysicalVolumes(G4VPhysicalVolume* vol1,
				      G4VPhysicalVolume* vol2);

	inline const G4VPhysicalVolume* GetVolume1() const;
	inline const G4VPhysicalVolume* GetVolume2() const;

	inline void SetVolume1(G4VPhysicalVolume* vol1);
	inline void SetVolume2(G4VPhysicalVolume* vol2);
          // These are potentially dangerous.

        static const G4RWTPtrOrderedVector<G4LogicalBorderSurface>* GetSurfaceTable();
        static size_t GetNumberOfBorderSurfaces();
	static void DumpInfo(); 
	  //   Methods dealing with the table of surfaces.

	inline size_t GetIndex() const;

        //////////////
        // Operators
        //////////////

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

	static G4RWTPtrOrderedVector<G4LogicalBorderSurface> theBorderSurfaceTable;
          // The static Table of Surfaces

	size_t theIndexInTable;
          // Index of surface in the surface table
};

typedef G4RWTPtrOrderedVector<G4LogicalBorderSurface> G4LogicalBorderSurfaceTable;

////////////////////
// Inline methods
////////////////////

#include "G4LogicalBorderSurface.icc"

#endif /* G4LogicalBorderSurface_h */
