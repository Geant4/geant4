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
// $Id: G4NURBShexahedron.hh,v 1.6 2001-07-11 10:01:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hexa hedron builder prototype
// OC 17 9 96

// Class Description:
// Hexa hedron builder prototype for NURBS.
// See documentation in graphics_reps/doc for details.
// Class Description - End:


#include "G4NURBS.hh"

#include "G4ThreeVector.hh"


#ifndef __C_G4NURBShexahedron__

#define __C_G4NURBShexahedron__ 1 

class	G4NURBShexahedron : public G4NURBS
	{
	
	// imagine the hexahedron is just a box, then
	// the eight corners must be given in the following order :
	//  DX  DY -DZ
	// -DX  DY -DZ
	// -DX -DY -DZ
	//  DX -DY -DZ
	//  DX  DY  DZ 
	// -DX  DY  DZ
	// -DX -DY  DZ
	//  DX -DY  DZ
 	// (ie, with Oz pointing to you, Ox on the right, Oy on the top:
	//  from the rear, from the upper right corner to the lower one
	// in anticlockwise sens, then the same for front side)

	public:
	  G4NURBShexahedron(const G4ThreeVector Corners [8]);
	  virtual G4Visible&  operator = (const G4Visible& right);
	  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
	  const char*	Whoami() const;
	};

#endif
// end of __C_G4NURBShexahedron__




