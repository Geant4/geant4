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
// $Id: G4OpticalSurface.cc,v 1.7 2001-10-17 07:59:54 gcosmo Exp $
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
		  			  theMaterialPropertiesTable(0)
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
        "  Surface type   = " << G4int(theType)   << G4endl <<
	"  Surface finish = " << G4int(theFinish) << G4endl <<
	"  Surface model  = " << G4int(theModel)  << G4endl;

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
