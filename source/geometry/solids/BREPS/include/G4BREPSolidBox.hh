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
// $Id: G4BREPSolidBox.hh,v 1.12 2010-10-20 09:14:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// Points have to be given in following way:
//  First for YZ plane with X fixed to negative value, then
//        for YZ plane with X xixed to positive value
//  Examle :     
//        const G4Point3D Pt1(-x,-y,-z);
//        const G4Point3D Pt2(-x,-y, z);
//        const G4Point3D Pt3(-x, y, z);
//        const G4Point3D Pt4(-x, y,-z);
//        const G4Point3D Pt5( x,-y,-z);
//        const G4Point3D Pt6( x,-y, z);
//        const G4Point3D Pt7( x, y, z);
//        const G4Point3D Pt8( x, y,-z);
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

  G4VSolid* Clone() const;
    // Returns a pointer of a dynamically allocated copy of the solid.

  std::ostream& StreamInfo(std::ostream& os) const;
    // Streams solid contents to output stream.

public:  // without description

  G4BREPSolidBox(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4BREPSolidBox(const G4BREPSolidBox& rhs);
  G4BREPSolidBox& operator=(const G4BREPSolidBox& rhs);
    // Copy constructor and assignment operator.

private:

  void InitializeBox();

private:

  G4RotationMatrix Rotation;
  G4Point3D constructorParams[8];

};

#endif
