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
// $Id$
//
// ----------------------------------------------------------------------
// Class G4CurvePoint
//
// Class Description:
//
// Class capable of storing both the parametric and the non-parametric
// representation of a point on a curve.
// The representation is evaluated lazily for efficiency.

// Author: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef included_G4CurvePoint
#define included_G4CurvePoint

#include "G4Curve.hh"
#include "G4Point3D.hh"

class G4CurvePoint 
{

public: // with description

  G4CurvePoint(G4Curve& c0);
    // Constructor, taking a curve as argument.

  virtual ~G4CurvePoint();
    // Empty destructor.

  G4CurvePoint(const G4CurvePoint& cp);
  G4CurvePoint& operator=(const G4CurvePoint& cp);
    // Copy constructor and assignment operator.

  inline void Init(G4Curve& c0);
    // Initialises a G4CurvePoint. Called by the constructor.

  inline G4Curve& GetCurve() const;
    // Returns the curve which the point belongs to.

  inline void Reset();
  inline void Reset(G4double u0);
  inline void Reset(const G4Point3D& p0);
  inline void Reset(G4double u0, const G4Point3D& p0);
    // Resets point's attributes.

  inline G4double GetPPoint();
  inline const G4Point3D& GetPoint();
    // Returns point as parameter or as point in space.

protected:
  
  G4CurvePoint();
    // Protected default constructor.

protected:

  // data

  G4Curve*  c;
  
  G4Point3D p;
  G4double  u;

  G4int notComputed;
  static const G4int pFlag;
  static const G4int uFlag;
  static const G4int allFlags;

};

#include "G4CurvePoint.icc"

#endif
