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
// $Id: G4BSplineSurfaceWithKnotsCreator.hh,v 1.5 2002-11-21 16:49:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BSplineSurfaceWithKnotsCreator
//
// Class description:
//
//

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------
#ifndef G4BSPLINESURFACEWITHKNOTSCREATOR_HH
#define G4BSPLINESURFACEWITHKNOTSCREATOR_HH

#include "G4BSplineSurfaceCreator.hh"

class G4BSplineSurfaceWithKnotsCreator: public G4BSplineSurfaceCreator
{
  public:

  // Constructor & destructor

    G4BSplineSurfaceWithKnotsCreator();
    ~G4BSplineSurfaceWithKnotsCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void*);
    const char* Name() const { return "B_Spline_Surface_With_Knots"; }
    static G4BSplineSurfaceWithKnotsCreator GetInstance();

  // Members

  private:

    static G4BSplineSurfaceWithKnotsCreator csc;
};

#endif
