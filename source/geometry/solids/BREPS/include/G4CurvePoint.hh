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
// $Id: G4CurvePoint.hh,v 1.5 2001-07-11 09:59:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
