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
// $Id: G4BREPSolidTorus.hh,v 1.7 2002-12-03 14:32:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidTorus
//
// Class description:
// 
//  Definition of a generic BREP torus:
//
//  G4BREPSolidTorus(const G4String& name,
//                   const G4ThreeVector& origin,
//                   const G4ThreeVector& axis,
//                   const G4ThreeVector& direction,
//                         G4double MinRadius,
//                         G4double MaxRadius)

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidTorus
#define __G4BREPSolidTorus

#include "G4BREPSolid.hh"

class G4BREPSolidTorus : public G4BREPSolid
{

 public:  // with description

  G4BREPSolidTorus(const G4String& name,
		   const G4ThreeVector& origin,
		   const G4ThreeVector& axis,
		   const G4ThreeVector& direction,
		   G4double MinRadius,
		   G4double MaxRadius);
    // Constructor defining the torus through its surfaces.

  ~G4BREPSolidTorus();
    // Empty destructor.

  virtual G4std::ostream& StreamInfo(G4std::ostream& os) const;
    // Streams solid contents to output stream.

 private:

  G4BREPSolidTorus(const G4BREPSolidTorus&);
  G4BREPSolidTorus& operator=(const G4BREPSolidTorus&);
    // Private copy constructor and assignment operator.

  struct G4BREPTorusParams
  {
    G4ThreeVector origin;
    G4ThreeVector axis;
    G4ThreeVector direction;
    G4double      MinRadius;
    G4double      MaxRadius;
  } constructorParams;
  
};

#endif
