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
// $Id: G4PlacedSolid.hh,v 1.8 2002-11-06 23:28:53 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PlacedSolid
//
// Class description:
//
// Class for a generic BREP solid placed in a 3D space.
// Attributes are: the solid, a traslation vector and a
// rotation matrix.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef G4PLACEDSOLID_HH
#define G4PLACEDSOLID_HH

#include "G4BREPSolid.hh"
#include "G4RotationMatrix.hh"

class G4PlacedSolid
{

public:  // with description

  G4PlacedSolid();
  G4PlacedSolid(G4BREPSolid*, G4Axis2Placement3D* =0);
  ~G4PlacedSolid();
    // Constructors & destructor

  inline G4bool operator==(const G4PlacedSolid& ps) const;
    // Equality operator.

  inline G4VSolid*         GetSolid() const;
  inline G4RotationMatrix* GetRotation() const;
  inline G4ThreeVector*    GetTranslation() const;
    // Accessors.

  G4std::ostream& StreamInfo(G4std::ostream& os) const;
    // Streams solid contents to output stream.

private:

  G4PlacedSolid(const G4PlacedSolid&);
  G4PlacedSolid& operator=(const G4PlacedSolid&);
    // Private copy constructor and assignment operator.

private:

  G4BREPSolid*      solid;
  G4RotationMatrix* solidRotation;
  G4ThreeVector*    solidTranslation;

};

#include "G4PlacedSolid.icc"

#endif
