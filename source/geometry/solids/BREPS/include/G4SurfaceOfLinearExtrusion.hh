// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SurfaceOfLinearExtrusion.hh,v 1.2 2000-08-28 08:57:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4SurfaceOfLinearExtrusion
//
// Class description:
// 
// Definition of a surface of linear extrusion.
// (not implemented yet).

// ----------------------------------------------------------------------
#ifndef included_G4SurfaceOfLinearExtrusion
#define included_G4SurfaceOfLinearExtrusion

#include "G4Surface.hh"

class G4SurfaceOfLinearExtrusion : public G4Surface
{

public:  // with description

  G4SurfaceOfLinearExtrusion();
  virtual ~G4SurfaceOfLinearExtrusion();
    // Constructor & destructor.

private:

  G4SurfaceOfLinearExtrusion(const G4SurfaceOfLinearExtrusion &);
  G4SurfaceOfLinearExtrusion& operator=(const G4SurfaceOfLinearExtrusion &);

};

#endif
