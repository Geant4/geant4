// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalSkinSurface.cc,v 1.2 1999-12-15 14:50:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// G4LogicalSkinSurface Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4LogicalSkinSurface.cc
// Description: A Logical Surface class for the surface  
//                surrounding a single logical volume.
// Version:     1.0
// Created:     1997-06-26
// Author:      John Apostolakis
// mail:        John.Apostolakis@cern.ch
// Modified:    1997-06-26  John Apostolakis
//
// CVS Id tag:  
////////////////////////////////////////////////////////////////////////

#include "G4LogicalSkinSurface.hh"
#include "G4ios.hh"
// #include "G4OpticalSurface.hh"

G4LogicalSkinSurfaceTable G4LogicalSkinSurface::theSurfaceTable;

/////////////////////////
// Class Implementation
/////////////////////////

  /////////////////
  // Constructors
  /////////////////

G4LogicalSkinSurface::G4LogicalSkinSurface(const G4String&   name,
					   G4LogicalVolume*  logicalVolume,
					   G4OpticalSurface* opticalSurface)
		: G4LogicalSurface(name, opticalSurface),
		  LogVolume(logicalVolume)
{
	// Store in the table of Surfaces
	theSurfaceTable.insert(this);
	theIndexInTable = theSurfaceTable.index(this);
}

G4LogicalSkinSurface::G4LogicalSkinSurface(const G4LogicalSkinSurface &right)
       : G4LogicalSurface(right.GetName(), right.GetOpticalSurface())
{
	*this = right;
}

G4LogicalSkinSurface::~G4LogicalSkinSurface(){}

  //////////////
  // Operators
  //////////////

const G4LogicalSkinSurface& G4LogicalSkinSurface::operator=(const G4LogicalSkinSurface &right)
{
	return right;
}

G4int G4LogicalSkinSurface::operator==(const G4LogicalSkinSurface &right) const
{
	return (this == (G4LogicalSkinSurface *) &right);
}

G4int G4LogicalSkinSurface::operator!=(const G4LogicalSkinSurface &right) const
{
	return (this != (G4LogicalSkinSurface *) &right);
}
  ////////////
  // Methods
  ////////////

G4LogicalSkinSurface* G4LogicalSkinSurface::GetSurface(const G4LogicalVolume* vol)
{
	for (int i=0; i<theSurfaceTable.length(); i++) {
		if(theSurfaceTable[i]->GetLogicalVolume() == vol)
			return theSurfaceTable[i];
	}
	return NULL;
}

void G4LogicalSkinSurface::DumpInfo() 
{

    // Dump info for known surfaces

    G4cout << "***** Surface Table : Nb of Surfaces = " << 
// G4LogicalSkinSurface::
      GetNumberOfSkinSurfaces() << " *****" << G4endl;

    for (int i=0; i<theSurfaceTable.length(); i++) {
      G4LogicalSkinSurface *pSkinSurface= theSurfaceTable[i];
      G4cout << theSurfaceTable[i]->GetName() << " : " << G4endl <<
	" Skin of logical volume " << pSkinSurface->GetLogicalVolume()->GetName  ()
	<< G4endl <<
	" Optical Surface Ptr  = " << (long) (pSkinSurface->GetOpticalSurface() )
	<< G4endl;

#ifdef PRINT_INFO
      //  DOES NOT COMPILE without including "G4OpticalSurface.hh"
    
      //  G4cout << pSkinSurface->GetOpticalSurface() << G4endl ;
      G4pticalSurface opticalSurface= pSkinSurface->GetOpticalSurface(); 
      G4cout << 
          "  Surface type   = " << opticalSurface->GetType()   << G4endl <<
          "  Surface finish = " << opticalSurface->GetFinish() << G4endl <<
	  "  Surface model  = " << opticalSurface->GetModel()  << G4endl;

/* 
        operator << ( G4OpticalSurface opticalSurface ) should exist 
	   and do something like:
          "  Surface type   = " << opticalSurface->GetType()   << G4endl <<
          "  Surface finish = " << opticalSurface->GetFinish() << G4endl <<
	  "  Surface model  = " << opticalSurface->GetModel()  << G4endl;
 */
#endif 

    }
    G4cout << G4endl;
}
