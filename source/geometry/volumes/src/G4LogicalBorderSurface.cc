// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalBorderSurface.cc,v 1.1 1999-01-07 16:08:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// G4LogicalBorderSurface Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4LogicalBorderSurface.cc
// Description: A Logical Surface class for surfaces defined by the
//               boundary of two physical volumes.
// Version:     1.0
// Created:     1997-06-26
// Author:      John Apostolakis
// mail:        John.Apostolakis@cern.ch
// Modified:    1997-06-26  John Apostolakis
//
// Id tag:      
////////////////////////////////////////////////////////////////////////

#include "G4LogicalBorderSurface.hh"


G4LogicalBorderSurfaceTable G4LogicalBorderSurface::theBorderSurfaceTable;

/////////////////////////
// Class Implementation
/////////////////////////

  /////////////////
  // Constructors
  /////////////////

G4LogicalBorderSurface::G4LogicalBorderSurface(const G4String& name,
					       G4VPhysicalVolume* vol1, 
					       G4VPhysicalVolume* vol2,
					       G4OpticalSurface* opticsSurface)
		: G4LogicalSurface(name, opticsSurface),
		  Volume1(vol1),
		  Volume2(vol2)
{
	// Store in the table of Surfaces
	theBorderSurfaceTable.insert(this);
	theIndexInTable = theBorderSurfaceTable.index(this);
}

G4LogicalBorderSurface::G4LogicalBorderSurface(const G4LogicalBorderSurface &right)
       : G4LogicalSurface(right.GetName(), right.GetOpticalSurface())
{
	*this = right;
}

G4LogicalBorderSurface::~G4LogicalBorderSurface(){}

  //////////////
  // Operators
  //////////////

const G4LogicalBorderSurface& G4LogicalBorderSurface::operator=(const G4LogicalBorderSurface &right)
{
	return right;
}

G4int G4LogicalBorderSurface::operator==(const G4LogicalBorderSurface &right) const
{
	return (this == (G4LogicalBorderSurface *) &right);
}

G4int G4LogicalBorderSurface::operator!=(const G4LogicalBorderSurface &right) const
{
	return (this != (G4LogicalBorderSurface *) &right);
}
  ////////////
  // Methods
  ////////////

G4LogicalBorderSurface* G4LogicalBorderSurface::GetSurface(const G4VPhysicalVolume* vol1,
				 const G4VPhysicalVolume* vol2)
{
	for (int i=0; i<theBorderSurfaceTable.length(); i++) {
		if(theBorderSurfaceTable[i]->GetVolume1() == vol1 &&
		   theBorderSurfaceTable[i]->GetVolume2() == vol2 )
			return theBorderSurfaceTable[i];
	}
	return NULL;
}

void G4LogicalBorderSurface::DumpInfo() // Class method (it is really const)
{

    // Dump info for known surfaces

    G4cout << "***** Surface Table : Nb of Surfaces = " << 
            GetNumberOfBorderSurfaces() << " *****" << endl;

    for (int i=0; i<theBorderSurfaceTable.length(); i++) {
      G4cout << theBorderSurfaceTable[i]->GetName() << " : " << endl <<
          "  Surface type   = " << theBorderSurfaceTable[i]->GetName() << endl;
#ifdef PRINT_INFO
          "  Surface type   = " << theBorderSurfaceTable[i]->GetOpticalSurface()->GetType()   << endl;
          "  Surface finish = " << theBorderSurfaceTable[i]->GetFinish() << endl <<
	  "  Surface model  = " << theBorderSurfaceTable[i]->GetModel()  << endl;
#endif 
    }
    G4cout << endl;
}
