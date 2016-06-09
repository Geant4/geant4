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
// Class G4BezierSurface
//
// Class description:
// 
// Definition of a rectangular trimmed surface.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __RECTRIMMEDSURFACE_H
#define __RECTRIMMEDSURFACE_H

#include "G4FCylindricalSurface.hh"

class G4RectangularTrimmedSurface : public G4Surface
{

public:  // with description

  G4RectangularTrimmedSurface();
  virtual ~G4RectangularTrimmedSurface();  
    // Constructor & destructor

  G4int Intersect(const G4Ray&);  
    // Counts the number of intersections of a the surface by a ray.
    // If no intersection is found, it sets distance = kInfinity and returns 0.

  void CalcBBox();
    // Computes the bounding-box.

  virtual const char* Name() const;
    // Returns the class type name.

private:

  G4RectangularTrimmedSurface(const G4RectangularTrimmedSurface&);
  G4RectangularTrimmedSurface& operator=(const G4RectangularTrimmedSurface&);

private:

  G4Surface* BasisSurface;

  G4double TrimU1,TrimU2;
  G4double TrimV1,TrimV2; 
 
  G4Point3D TrimPointU1, TrimPointU2;
  G4Point3D TrimPointV1, TrimPointV2;  
  
};

#endif
