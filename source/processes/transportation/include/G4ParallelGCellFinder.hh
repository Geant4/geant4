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
// $Id: G4ParallelGCellFinder.hh,v 1.1 2002-10-16 16:27:41 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelGCellFinder
//
// Class description:
//
// Find a G4GeometryCell in the parallel geometry.
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4ParallelGCellFinder_hh
#define G4ParallelGCellFinder_hh G4ParallelGCellFinder_hh

#include "globals.hh"

#include "G4VGCellFinder.hh"

class G4VParallelStepper;
class  G4ParallelGCellFinder : public G4VGCellFinder 
{

public:  // with description

  explicit G4ParallelGCellFinder(const G4VParallelStepper &astepper);
  virtual  ~G4ParallelGCellFinder();

  virtual G4GeometryCell 
  GetPreGeometryCell(const G4Step &aStep) const;
  virtual G4GeometryCell 
  GetPostGeometryCell(const G4Step &aStep) const;

private:  
  const G4VParallelStepper &fPStepper;
  
};

#endif
