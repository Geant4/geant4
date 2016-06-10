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
// $Id: G4VGCellFinder.hh 66356 2012-12-18 09:02:32Z gcosmo $
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
