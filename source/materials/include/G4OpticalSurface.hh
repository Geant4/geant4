// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpticalSurface.hh,v 1.2 1999-04-14 12:49:00 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// G4OpticalSurface Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpticalSurface.hh
// Description: A optical surface class for use in G4OpBoundaryProcess
// Version:     2.0
// Created:     1997-06-26
// Author:      Peter Gumplinger
// mail:        gum@triumf.ca
//
// Cvs version: 
////////////////////////////////////////////////////////////////////////

#ifndef G4OpticalSurface_h
#define G4OpticalSurface_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"

enum G4OpticalSurfaceFinish
{
   polished,                    // smooth perfectly polished surface
   polishedfrontpainted,        // smooth top-layer (front) paint
   polishedbackpainted,         // same is 'polished' but with a back-paint
   ground,                      // rough surface
   groundfrontpainted,          // rough top-layer (front) paint
   groundbackpainted            // same as 'ground' but with a back-paint
};

enum G4OpticalSurfaceType
{
   dielectric_metal,            // dielectric-metal interface
   dielectric_dielectric        // dielectric-dielectric interface
};

enum G4OpticalSurfaceModel
{
   glisur,                      // original GEANT3 model
   unified                      // UNIFIED model
};

class G4MaterialPropertiesTable;

/////////////////////
// Class Definition
/////////////////////

class G4OpticalSurface
{

public:

        //////////////
        // Operators
        //////////////

        G4OpticalSurface(const G4OpticalSurface &right);
	const G4OpticalSurface & operator=(const G4OpticalSurface &right);

	G4int operator==(const G4OpticalSurface &right) const;
	G4int operator!=(const G4OpticalSurface &right) const;

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

	G4OpticalSurface(const G4String& name,
			 G4OpticalSurfaceModel model = glisur,
			 G4OpticalSurfaceFinish finish = polished,
			 G4OpticalSurfaceType type = dielectric_dielectric,
			 G4double value = 1.0);

	~G4OpticalSurface();

	////////////
	// Methods
        ////////////

public:

	// public methods

	G4String GetName() const { return theName; };
	void     SetName(const G4String& name){theName = name;};

        // set/get the surface type

        G4OpticalSurfaceType GetType() const {return theType;};
        void         SetType(const G4OpticalSurfaceType type){theType = type;};

        // set/get the surface finish

        G4OpticalSurfaceFinish GetFinish() const {return theFinish;};
        void         SetFinish(const G4OpticalSurfaceFinish finish)
						 {theFinish = finish;};

        // set/get the surface model to be followed (glisur || unified)

        G4OpticalSurfaceModel GetModel() const {return theModel;};
        void           SetModel(const G4OpticalSurfaceModel model)
						   {theModel = model;};

	G4double GetSigmaAlpha() const {return sigma_alpha;};
	void     SetSigmaAlpha(const G4double s_a)
				        {sigma_alpha = s_a;};

	G4double GetPolish() const {return polish;};
	void     SetPolish(const G4double plsh) {polish=plsh;};

	void SetMaterialPropertiesTable(G4MaterialPropertiesTable *anMPT)
				    { theMaterialPropertiesTable = anMPT;};

	G4MaterialPropertiesTable* GetMaterialPropertiesTable() const
				       { return theMaterialPropertiesTable;};

	void DumpInfo() const;

private:

// ------------------
// Basic data members ( To define an optical surface)
// ------------------

	G4String theName;			// Surface name

        G4OpticalSurfaceModel theModel;		// Surface model
        G4OpticalSurfaceFinish theFinish;	// Surface finish
	G4OpticalSurfaceType theType;		// Surface type

	G4double sigma_alpha;		// The sigma of micro-facet polar angle
	G4double polish;		// Polish parameter in glisur model

	G4MaterialPropertiesTable* theMaterialPropertiesTable;

};

////////////////////
// Inline methods
////////////////////

#endif /* G4OpticalSurface_h */
