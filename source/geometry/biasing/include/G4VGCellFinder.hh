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
// $Id: G4VGCellFinder.hh,v 1.1 2002-10-16 14:30:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VGCellFinder
//
// Class description:
//
// An interface of an GCellFinder. The G4GeometryCell is obtained
// in the parallel geometry from G4VParallelStepper and in the
// mass geometry form G4Step. This interface allows to implement
// the two different ways to get a G4GeometryCell used by a process.
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VGCellFinder_hh
#define G4VGCellFinder_hh G4VGCellFinder_hh

#include "globals.hh"
#include "G4GeometryCell.hh"

class G4Step;

class  G4VGCellFinder
{

public:  // with description

  G4VGCellFinder();
  virtual  ~G4VGCellFinder();

  virtual G4GeometryCell 
  GetPreGeometryCell(const G4Step &aStep) const=0;
  virtual G4GeometryCell 
  GetPostGeometryCell(const G4Step &aStep) const=0;
  
};

#endif
