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
// $Id: G4SurfaceProperty.hh,v 1.1 2003-11-28 00:28:36 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// G4SurfaceProperty Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4SurfaceProperty.hh
// Description: A base class for for descriping surface property such
//              as G4OpticalSurface, G4FirsovSurface, G4X-raySurface
// Version:     1.0
// Created:     13-10-2003
// Author:      Fan Lei
// Updated:     27-11-2003: add method and class descriptors
// mail:        fan.lei@cern.ch
//
// Cvs version: 
////////////////////////////////////////////////////////////////////////

#ifndef G4SurfaceProperty_h
#define G4SurfaceProperty_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"

// Class Description:
// A base class for descriping a surface property. Derived classes are G4Opticalsurface,
// G4Firovsurface, etc.      
// Contains the enumerations: G4SurfaceType
// Class Description - End:

enum G4SurfaceType
{
   dielectric_metal,            // dielectric-metal interface
   dielectric_dielectric,       // dielectric-dielectric interface
   firsov,                      // for Firsov Process
   x_ray                        // for x-ray mirror process
};

/////////////////////
// Class Definition
/////////////////////

class G4SurfaceProperty
{

public: // Without description

        //////////////
        // Operators
        //////////////

  //        G4SurfaceProperty(const G4SurfaceProperty &right);
  //	const G4SurfaceProperty & operator=(const G4SurfaceProperty &right);

  //	G4int operator==(const G4SurfaceProperty &right) const;
  //	G4int operator!=(const G4SurfaceProperty &right) const;

public: // With description

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        G4SurfaceProperty(const G4String& name,
			  G4SurfaceType type = x_ray)
	  :theName(name), theType(type)
        {
        }
        // Constructor of a X-ray optical surface object.

public: // Without description

	~G4SurfaceProperty(){};

	////////////
	// Methods
        ////////////

	// public methods

public: // With description

	G4String GetName() const { return theName; };
        // Returns the surface name.
	void     SetName(const G4String& name){theName = name;};
        // Sets the surface name.

        G4SurfaceType GetType() const {return theType;};
        // Returns the surface type.
        void         SetType(const G4SurfaceType type){theType = type;};
        // Sets the surface type.        

protected:

// ------------------
// Basic data members ( To define surface property)
// ------------------

	G4String theName;			// Surface name

	G4SurfaceType theType;		// Surface type

};

////////////////////
// Inline methods
////////////////////

#endif /* G4SurfaceProperty_h */
