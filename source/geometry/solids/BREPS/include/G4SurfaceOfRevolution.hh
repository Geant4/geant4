// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SurfaceOfRevolution.hh,v 1.3 2000-11-08 14:22:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4SurfaceOfRevolution
//
// Class description:
// 
// Definition of a surface of revolution.
// (not implemented yet).

// ----------------------------------------------------------------------
#ifndef included_G4SurfaceOfRevolution
#define included_G4SurfaceOfRevolution

#include "G4Surface.hh"

class G4SurfaceOfRevolution : public G4Surface
{

public:  // with description

  G4SurfaceOfRevolution();
  virtual ~G4SurfaceOfRevolution();
    // Constructor & destructor.

private:

  G4SurfaceOfRevolution(const G4SurfaceOfRevolution &);
  G4SurfaceOfRevolution& operator=(const G4SurfaceOfRevolution &);
    // Private copy constructor and assignment operator.
};

#endif
