// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpticalSurface.cc,v 1.4 1999-12-15 14:50:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// Optical Surface Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpticalSurface.cc
// Description: An optical surface class for use in G4OpBoundaryProcess
// Version:     2.0
// Created:     1997-06-26
// Author:      Peter Gumplinger
// mail:        gum@triumf.ca
//
// Cvs version: 
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4OpticalSurface.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

const G4OpticalSurface& 
      G4OpticalSurface::operator=(const G4OpticalSurface& right)
{
  if (this != &right)
    {
      theName                    = right.theName;
      theModel                   = right.theModel;
      theFinish                  = right.theFinish;
      theType                    = right.theType;
      sigma_alpha                = right.sigma_alpha;
      polish                     = right.polish;
      theMaterialPropertiesTable = right.theMaterialPropertiesTable;
     } 
  return *this;
}

        /////////////////
        // Constructors
        /////////////////

G4OpticalSurface::G4OpticalSurface(const G4String& name,
				   G4OpticalSurfaceModel model,
				   G4OpticalSurfaceFinish finish,
				   G4OpticalSurfaceType type,
				   G4double value)
					: theName(name),
		  			  theModel(model),
		  			  theFinish(finish),
		  			  theType(type),
		  			  theMaterialPropertiesTable(NULL)
{
	if (model == glisur ){
		polish = value;
		sigma_alpha = 0.0;
	}
	else if ( model == unified ) {
		sigma_alpha = value;
		polish = 0.0;
	}
	else {
		G4Exception("G4OpticalSurface::G4OpticalSurface ==> " 
			    "Constructor called with INVALID model.");
	}
}

G4OpticalSurface::G4OpticalSurface(const G4OpticalSurface &right)
{
	*this = right;
}

G4OpticalSurface::~G4OpticalSurface(){}

G4int G4OpticalSurface::operator==(const G4OpticalSurface &right) const
{
	return (this == (G4OpticalSurface *) &right);
}

G4int G4OpticalSurface::operator!=(const G4OpticalSurface &right) const
{
	return (this != (G4OpticalSurface *) &right);
}
        ////////////
        // Methods
        ////////////

void G4OpticalSurface::DumpInfo() const 
{

	// Dump info for surface

	G4cout << 
        "  Surface type   = " << theType   << G4endl <<
	"  Surface finish = " << theFinish << G4endl <<
	"  Surface model  = " << theModel  << G4endl;

	G4cout << G4endl;

	G4cout << "  Surface parameter " << G4endl;
	G4cout << "  ----------------- " << G4endl;
	if (theModel == glisur ){
		G4cout << polish      << G4endl;
	}
	else {
		G4cout << sigma_alpha << G4endl;
	}
	G4cout << G4endl;
}
