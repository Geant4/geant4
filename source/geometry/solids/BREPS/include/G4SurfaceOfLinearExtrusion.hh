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
// $Id: G4SurfaceOfLinearExtrusion.hh,v 1.4 2001-07-11 09:59:38 gunter Exp $
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
    // Private copy constructor and assignement operator.

};

#endif
