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
// $Id: G4PointRat.hh,v 1.10 2006-06-29 18:40:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PointRat
//
// Class description:
//
// A G4PointRat object is composed of:
// - a point in 3D space (G4Point3D)
// - a scale factor (set to 1 by default)

// Authors: J.Sulkimo, P.Urban.
// Revisions by: A.Floquet, L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4POINT_RAT
#define __G4POINT_RAT

#include "geomdefs.hh"              // For kInfinity

#include "G4Point3D.hh"
#include "G4Plane3D.hh"

#define SQRT_SMALL_FASTF 1.0e-18
#define SMALL            SQRT_SMALL_FASTF       
#define ROW              0
#define COL              1

const G4Point3D PINFINITY(kInfinity, kInfinity, kInfinity);

class G4PointRat
{

public:  // with description

  G4PointRat();
  ~G4PointRat();
    // Constructor & destructor.

  G4PointRat(const G4Point3D&);
    // Copy constructor.

  G4PointRat& operator=(const G4Point3D&);
  G4PointRat& operator=(const G4PointRat&);
    // Overloaded assignment operators.

  inline G4double x() const;
  inline void setX (G4double Value);
  inline G4double y() const;
  inline void setY (G4double Value);
  inline G4double z() const;
  inline void setZ (G4double Value);
  inline G4double w() const;
  inline void setW(G4double Value);
  inline G4Point3D pt() const;
    // Get/Set methods for attributes.

  inline G4double PlaneDistance(const G4Plane3D& Pl) const;
    // Returns distance from a given plane.

public:  // without description

  G4int GetType(void) const; // This function should be removed
                             // if calls to this are also removed

private:

   G4Point3D pt3d;
   G4double  s;
};

#include "G4PointRat.icc"

#endif
