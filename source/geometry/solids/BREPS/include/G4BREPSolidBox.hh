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
// $Id: G4BREPSolidBox.hh,v 1.8 2005/11/09 15:01:25 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidBox
//
// Class description:
//
//  Definition of a generic BREP box:
//
//  G4BREPSolidBox(const G4String& name,
//                 const G4Point3D& Pt1, const G4Point3D& Pt2, 
//                 const G4Point3D& Pt3, const G4Point3D& Pt4,
//                 const G4Point3D& Pt5, const G4Point3D& Pt6,
//                 const G4Point3D& Pt7, const G4Point3D& Pt8)      

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidBOX
#define __G4BREPSolidBOX

#include "G4BREPSolid.hh"
#include "G4RotationMatrix.hh"

class G4BREPSolidBox : public G4BREPSolid
{

public: // with description

  G4BREPSolidBox(const G4String&,  const G4Point3D&, const G4Point3D&, 
		 const G4Point3D&, const G4Point3D&, const G4Point3D&, 
		 const G4Point3D&, const G4Point3D&, const G4Point3D& );       
    // Constructor defining the box through its planes.

  ~G4BREPSolidBox();
    // Empty destructor.

  EInside Inside(register const G4ThreeVector& Pt) const;
    // Determines if the point Pt is inside, outside or on the surface
    // of the solid.

  virtual std::ostream& StreamInfo(std::ostream& os) const;
    // Streams solid contents to output stream.

public:  // without description

  G4BREPSolidBox(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

private:

  G4BREPSolidBox(const G4BREPSolidBox&);
  G4BREPSolidBox& operator=(const G4BREPSolidBox&);
    // Private copy constructor and assignment operator.

private:

  G4RotationMatrix Rotation;
  G4Point3D constructorParams[8];

};

#endif
