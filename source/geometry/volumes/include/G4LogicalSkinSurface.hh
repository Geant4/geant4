// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalSkinSurface.hh,v 1.4 2000-04-25 16:15:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// class G4LogicalSkinSurface
////////////////////////////////////////////////////////////////////////
//
// Class description:
//
// A Logical Surface class for the surface surrounding a single logical
// volume.

// History:
// -------
// Created:     1997-06-16
// Author:      John Apostolakis
// mail:        John.Apostolakis@cern.ch
// Modified:    1997-06-16  John Apostolakis
//
////////////////////////////////////////////////////////////////////////

#ifndef G4LogicalSkinSurface_h
#define G4LogicalSkinSurface_h 1

/////////////
// Includes
/////////////

#include "G4LogicalSurface.hh"
#include "G4LogicalVolume.hh"

// G4RWTPtrOrderedVector
#include "g4rw/tpordvec.h"


/////////////////////
// Class Definition
/////////////////////


class G4LogicalSkinSurface : public G4LogicalSurface 
{

  public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        // Is the name meaningful for the logical skin surface ?

        G4LogicalSkinSurface(const G4String& name, G4LogicalVolume* vol,
			     G4OpticalSurface* opticalSurface);

	~G4LogicalSkinSurface();

	////////////
	// Methods
        ////////////
  public:
  
        static G4LogicalSkinSurface* GetSurface(const G4LogicalVolume* vol);

	G4LogicalVolume* GetLogicalVolume() const;
	void		 SetLogicalVolume(G4LogicalVolume* vol);

        static size_t GetNumberOfSkinSurfaces();
        static void DumpInfo(); // const 
	  // Methods dealing with the table of surfaces.

#if THESE_ARE_NEEDED
	size_t GetIndex() const;
#endif

        //////////////
        // Operators
        //////////////

	G4int operator==(const G4LogicalSkinSurface &right) const;
	G4int operator!=(const G4LogicalSkinSurface &right) const;

private:

        G4LogicalSkinSurface(const G4LogicalSkinSurface &right);
	const G4LogicalSkinSurface & operator=(const G4LogicalSkinSurface &right);
	  // Assignment and copying must be denied.

	// ------------------
	// Basic data members ( To define a 'logical' surface)
	// ------------------

private:
	G4LogicalVolume* LogVolume;
	  // Logical Volume pointer on side 1

  	static G4RWTPtrOrderedVector<G4LogicalSkinSurface> theSurfaceTable;
	  // The static Table of Surfaces

	size_t theIndexInTable;
	  // Index of surface in the surface table

};

typedef G4RWTPtrOrderedVector<G4LogicalSkinSurface> G4LogicalSkinSurfaceTable;

////////////////////
// Inline methods
////////////////////

#include "G4LogicalSkinSurface.icc"

#endif /* G4LogicalSkinSurface_h */
