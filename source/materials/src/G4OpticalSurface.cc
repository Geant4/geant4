//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4OpticalSurface.cc,v 1.10 2006/06/29 19:13:08 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "globals.hh"
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
      theName                    = right.GetName();
      theModel                   = right.theModel;
      theFinish                  = right.theFinish;
      theType                    = right.GetType();
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
				   G4SurfaceType type,
				   G4double value)
                                   : G4SurfaceProperty(name,type),
				     theModel(model),
				     theFinish(finish),
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
  : G4SurfaceProperty(right.GetName())
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
