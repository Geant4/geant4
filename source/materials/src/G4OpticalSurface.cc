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
// $Id: G4OpticalSurface.cc,v 1.17 2010-10-25 15:16:02 vnivanch Exp $
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
      AngularDistribution        = right.AngularDistribution;
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
        else if ( model == LUT ) {
                sigma_alpha = value;
                polish = 0.0;
        }
	else {
		G4Exception("G4OpticalSurface::G4OpticalSurface ==> " 
			    "Constructor called with INVALID model.");
	}

        AngularDistribution = NULL;

        if (type == dielectric_LUT) {
           AngularDistribution =
                       new G4float[incidentIndexMax*thetaIndexMax*phiIndexMax];
           ReadFile();
        }
}

G4OpticalSurface::~G4OpticalSurface()
{
        if (AngularDistribution) delete AngularDistribution;
}

G4OpticalSurface::G4OpticalSurface(const G4OpticalSurface &right)
  : G4SurfaceProperty(right.GetName())
{
	*this = right;
}

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
        else if (theModel == LUT ){
                G4cout << sigma_alpha << G4endl;
        }
	else {
		G4cout << sigma_alpha << G4endl;
	}
	G4cout << G4endl;
}

void G4OpticalSurface::SetType(const G4SurfaceType& type)
{
  theType = type;
  if (type == dielectric_LUT) {
     if (!AngularDistribution) AngularDistribution =
                       new G4float[incidentIndexMax*thetaIndexMax*phiIndexMax];
     ReadFile();
  }
}

void G4OpticalSurface::SetFinish(const G4OpticalSurfaceFinish finish)
{
  theFinish = finish;
  if (theType == dielectric_LUT) {
     if (!AngularDistribution) AngularDistribution =
                       new G4float[incidentIndexMax*thetaIndexMax*phiIndexMax];
     ReadFile();
  }
}

void G4OpticalSurface::ReadFile()
{
  G4String readFileName = " ";

  if (theFinish == polishedlumirrorglue) {
     readFileName = "PolishedLumirrorGlue.dat";
  }
  else if (theFinish == polishedlumirrorair) {
     readFileName = "PolishedLumirror.dat";
  }
  else if (theFinish == polishedteflonair) {
     readFileName = "PolishedTeflon.dat";
  }
  else if (theFinish == polishedtioair) {
     readFileName = "PolishedTiO.dat";
  }
  else if (theFinish == polishedtyvekair) {
     readFileName = "PolishedTyvek.dat";
  }
  else if (theFinish == polishedvm2000glue) {
     readFileName = "PolishedVM2000Glue.dat";
  }
  else if (theFinish == polishedvm2000air) {
     readFileName = "PolishedVM2000.dat";
  }
  else if (theFinish == etchedlumirrorglue) {
     readFileName = "EtchedLumirrorGlue.dat";
  }
  else if (theFinish == etchedlumirrorair) {
     readFileName = "EtchedLumirror.dat";
  }
  else if (theFinish == etchedteflonair) {
     readFileName = "EtchedTeflon.dat";
  }
  else if (theFinish == etchedtioair) {
     readFileName = "EtchedTiO.dat";
  }
  else if (theFinish == etchedtyvekair) {
     readFileName = "EtchedTyvek.dat";
  }
  else if (theFinish == etchedvm2000glue) {
     readFileName = "EtchedVM2000Glue.dat";
  }
  else if (theFinish == etchedvm2000air) {
     readFileName = "EtchedVM2000.dat";
  }
  else if (theFinish == groundlumirrorglue) {
     readFileName = "GroundLumirrorGlue.dat";
  }
  else if (theFinish == groundlumirrorair) {
     readFileName = "GroundLumirror.dat";
  }
  else if (theFinish == groundteflonair) {
     readFileName = "GroundTeflon.dat";
  }
  else if (theFinish == groundtioair) {
     readFileName = "GroundTiO.dat";
  }
  else if (theFinish == groundtyvekair) {
     readFileName = "GroundTyvek.dat";
  }
  else if (theFinish == groundvm2000glue) {
     readFileName = "GroundVM2000Glue.dat";
  }
  else if (theFinish == groundvm2000air) {
     readFileName = "GroundVM2000.dat";
  }

  if (readFileName == " ") return;

  char* path = getenv("G4REALSURFACEDATA");
  if (!path) {
     G4String excep =
        "G4OpBoundaryProcess - G4REALSURFACEDATA environment variable not set";
     G4Exception(excep);
     return;
  }
  G4String pathString(path);

  readFileName = pathString + "/" + readFileName;

  // Open LUT with Material and Integer Angle
  FILE* readFileHandle;

  readFileHandle = fopen(readFileName,"r");

  if (readFileHandle) {
     G4int idxmax = incidentIndexMax*thetaIndexMax*phiIndexMax;
     for (G4int i=0;i<idxmax;i++) {
       fscanf(readFileHandle,"%6f", &AngularDistribution[i]);
     }
     G4cout << "LUT - data file: " << readFileName << " read in! " << G4endl;
  }
  else {
     G4String excep = "LUT - data file: " + readFileName + " not found";
     G4Exception(excep);
     return;
  }
  fclose(readFileHandle);
}
