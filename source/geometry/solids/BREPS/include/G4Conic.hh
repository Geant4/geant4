// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Conic.hh,v 1.4 2000-08-28 15:00:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Conic
//
// Class description:
// 
// Definition of a generic conical curve.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __CONIC_H
#define __CONIC_H 

#include "G4Curve.hh"
#include "G4Axis2Placement3D.hh"

class G4Conic : public G4Curve
{

public:  // with description
 
  G4Conic();
  virtual ~G4Conic();
    // Constructor & destructor.

  inline const G4Axis2Placement3D* GetPosition() const;
  inline G4double GetPShift() const;
  inline void SetPShift(G4double pShift0); 
    // Get/Set geometric data.

public:  // without description

  //inline G4Placement GetPosition() {return Position;}
  //virtual const char *Name(){return "G4ConicalCurve";}

protected:

  //void ProjectCurve(const G4Plane&, const G4Plane&);
  //G4int HitPartOfCurve(G4double, G4double, const G4Point2d&);
  //G4Placement Position;

  G4Axis2Placement3D position;
    // Geometric data.

private:

  G4double pShift;
    // pShift must be added/subtracted from the parameter.
    // No STEP I/O if not 0!!!
    // Set by Project members.

};

#include "G4Conic.icc"

#endif
