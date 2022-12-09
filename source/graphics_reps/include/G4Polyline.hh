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
//
// 
// John Allison  July 1995

// Class Description:
// A set of line segments defined with a set of vertices.
// G4Polyline is used for visualizing trajectories, steps, coordinate axes,
// etc.
// Class Description - End:


#ifndef G4POLYLINE_HH
#define G4POLYLINE_HH

#include "G4Visible.hh"
#include "G4Point3DList.hh"

#include "G4Transform3D.hh"

class G4Polyline: public G4Visible, public G4Point3DList {

  friend std::ostream& operator << (std::ostream& os, const G4Polyline&);

public: // With description

  using iterator = G4Point3DList::iterator;

  G4Polyline ();
  ~G4Polyline () override;
  G4Polyline& transform (const G4Transform3D&);
};

#endif
