// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalSkinSurface.hh,v 1.1 1999-01-07 16:08:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// G4LogicalSkinSurface Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4LogicalSkinSurface.hh
// Description: A Logical Surface class for the surface  
//                surrounding a single logical volume.
// Version:     1.0
// Created:     1997-06-16
// Author:      John Apostolakis
// mail:        John.Apostolakis@cern.ch
// Modified:    1997-06-16  John Apostolakis
//
// Id tag:      
////////////////////////////////////////////////////////////////////////

#ifndef G4LogicalSkinSurface_h
#define G4LogicalSkinSurface_h 1

/////////////
// Includes
/////////////

#include "G4LogicalSurface.hh"
#include "G4LogicalVolume.hh"

// RWTPtrOrderedVector
#include <rw/tpordvec.h>


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
	// public methods

        static G4LogicalSkinSurface* GetSurface(const G4LogicalVolume* vol);

	G4LogicalVolume* GetLogicalVolume() const {return LogVolume;};
	void		 SetLogicalVolume(G4LogicalVolume* vol)
						{LogVolume = vol;};


	//   Methods dealing with the table of surfaces.
	//
        static size_t GetNumberOfSkinSurfaces()
                                      { return theSurfaceTable.length(); };
        static void DumpInfo(); // const 

#if THESE_ARE_NEEDED
	// static const RWTPtrOrderedVector<G4LogicalSkinSurface>* 
	//          GetSurfaceTable()
	//                              { return &theSurfaceTable; };

	size_t GetIndex() const { return theIndexInTable; };
#endif

        //////////////
        // Operators
        //////////////
public:
	G4int operator==(const G4LogicalSkinSurface &right) const;
	G4int operator!=(const G4LogicalSkinSurface &right) const;

private:
	// Assignment and copying must be denied.
        G4LogicalSkinSurface(const G4LogicalSkinSurface &right);
	const G4LogicalSkinSurface & operator=(const G4LogicalSkinSurface &right);

	// ------------------
	// Basic data members ( To define a 'logical' surface)
	// ------------------
private:
	G4LogicalVolume* LogVolume;	// Logical Volume pointer on side 1

//	The static Table of Surfaces
  	static RWTPtrOrderedVector<G4LogicalSkinSurface> theSurfaceTable;
//	static G4LogicalSkinSurfaceTable theSurfaceTable;

	size_t theIndexInTable;		// Index of surface in the surface table

};

typedef RWTPtrOrderedVector<G4LogicalSkinSurface> G4LogicalSkinSurfaceTable;

////////////////////
// Inline methods
////////////////////

#endif /* G4LogicalSkinSurface_h */
