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
// $Id: G4LogicalSkinSurface.hh,v 1.7 2001-07-11 10:00:28 gunter Exp $
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

#include "g4std/vector"

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

	inline const G4LogicalVolume* GetLogicalVolume() const;
	inline void  SetLogicalVolume(G4LogicalVolume* vol);

        static size_t GetNumberOfSkinSurfaces();
        static void DumpInfo(); // const 
	  // Methods dealing with the table of surfaces.

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

  	static G4std::vector<G4LogicalSkinSurface*> theSurfaceTable;
	  // The static Table of Surfaces

};

typedef G4std::vector<G4LogicalSkinSurface*> G4LogicalSkinSurfaceTable;

////////////////////
// Inline methods
////////////////////

#include "G4LogicalSkinSurface.icc"

#endif /* G4LogicalSkinSurface_h */
