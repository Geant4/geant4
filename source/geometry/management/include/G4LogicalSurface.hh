// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalSurface.hh,v 1.4 2000-11-01 15:39:32 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// Class G4LogicalSurface
////////////////////////////////////////////////////////////////////////
//
// Class description:
//
// An abstraction of a geometrical surface, it is an abstract 
// base class for different implementations of surfaces.
// Its primary function is to hold pointers to objects that describe the
// surface's physical properties. For example it holds a pointer to a
// surface's optical properties, and because of this it is used in processes
// like G4OpBoundaryProcess. 
//
// Methods:
//   G4OpticalSurface*  GetOpticalSurface() const
//   void     SetOpticalSurface(G4OpticalSurface*)
//
//   G4String GetName() const
//   void     SetName(const G4String&)
//
//   G4TransitionRadiationSurface*  GetTransitionRadiationSurface() const
//   void SetTransitionRadiationSurface(G4TransitionRadiationSurface*)
//
// Data members:
//   G4String                       theName
//   G4OpticalSurface*              theOpticalSurface
//   G4TransitionRadiationSurface*  theTransRadSurface

// Version:     1.0
// Created:     1997, June, 4th to 17th
// Author:      John Apostolakis, (with help of Peter Gumplinger)
// mail:        japost@mail.cern.ch
// Modified:    1997, June 26th  John Apostolakis
//
// Id tag:      
////////////////////////////////////////////////////////////////////////

#ifndef G4LogicalSurface_h
#define G4LogicalSurface_h 1

/////////////
// Includes
/////////////

#include "globals.hh"

class G4OpticalSurface;
class G4TransitionRadiationSurface;

/////////////////////
// Class Definition
/////////////////////

class G4LogicalSurface
{
	////////////
	// Methods
        ////////////
 public:
   G4OpticalSurface*  GetOpticalSurface() const;
   void     SetOpticalSurface(G4OpticalSurface* ptrOpticalSurface);

   G4String GetName() const;
   void     SetName(const G4String& name);

   G4TransitionRadiationSurface*  GetTransitionRadiationSurface() const;
   void SetTransitionRadiationSurface(G4TransitionRadiationSurface* transRadSurf);

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

 protected:

        // There should be no instances of this class

   G4LogicalSurface(const G4String&         name,
		    G4OpticalSurface* opticalSurface); 
        // Is the name more meaningful for the properties or the logical
        // surface ?  

 public:
   virtual ~G4LogicalSurface();

 private:
   G4LogicalSurface(const G4LogicalSurface &right); // Copying restricted

        //////////////
        // Operators
        //////////////
 public:
   G4int operator==(const G4LogicalSurface &right) const;
   G4int operator!=(const G4LogicalSurface &right) const;

 private:
   const G4LogicalSurface& operator=(const G4LogicalSurface& right);

	// ------------------
	// Basic data members ( To define a 'logical' surface)
	// ------------------

 private:
   G4String theName;		// Surface name

   G4OpticalSurface*              theOpticalSurface;
   G4TransitionRadiationSurface*  theTransRadSurface;
};

////////////////////
// Inline methods
////////////////////

#include "G4LogicalSurface.icc"

#endif /* G4LogicalSurface_h */
