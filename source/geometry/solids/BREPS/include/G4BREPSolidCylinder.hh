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
// $Id: G4BREPSolidCylinder.hh,v 1.6 2002-11-06 23:28:47 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidCylinder
//
// Class description:
//
//  Definition of a generic BREP cylinder:
//
//  G4BREPSolidCylinder(const G4String& name,
//                      const G4ThreeVector& origin,
//                      const G4ThreeVector& axis,
//                      const G4ThreeVector& direction,		      
//                            G4double radius,
//                            G4double length)

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidCylinder
#define __G4BREPSolidCylinder

#include "G4BREPSolid.hh"

class G4BREPSolidCylinder : public G4BREPSolid
{

 public: // with description

  G4BREPSolidCylinder(const G4String& name,
		      const G4ThreeVector& origin,
		      const G4ThreeVector& axis,
		      const G4ThreeVector& direction,		      
		      G4double radius,
		      G4double length);
    // Constructor defining the cylinder through planes
    // and circular surfaces.

  ~G4BREPSolidCylinder();
    // Empty destructor.

  virtual G4std::ostream& StreamInfo(G4std::ostream& os) const;
    // Streams solid contents to output stream.

 private:

    G4BREPSolidCylinder(const G4BREPSolidCylinder&);
    G4BREPSolidCylinder& operator=(const G4BREPSolidCylinder&);
      // Private copy constructor and assignment operator.

  struct {
    G4ThreeVector origin;
    G4ThreeVector axis;
    G4ThreeVector direction;
    G4double      length;
    G4double      radius;
  } constructorParams;
};

#endif
