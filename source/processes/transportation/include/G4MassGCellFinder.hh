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
// $Id: G4MassGCellFinder.hh,v 1.1 2002-10-16 16:27:41 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelGCellFinder
//
// Class description:
//
// Find a G4GeometryCell in the mass geometry.
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4MassGCellFinder_hh
#define G4MassGCellFinder_hh G4MassGCellFinder_hh

#include "globals.hh"

#include "G4VGCellFinder.hh"

class  G4MassGCellFinder : public G4VGCellFinder 
{

public:  // with description

  G4MassGCellFinder();
  virtual  ~G4MassGCellFinder();

  virtual G4GeometryCell 
  GetPreGeometryCell(const G4Step &aStep) const;
  virtual G4GeometryCell 
  GetPostGeometryCell(const G4Step &aStep) const;

};

#endif
